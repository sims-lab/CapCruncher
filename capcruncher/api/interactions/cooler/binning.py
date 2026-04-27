from __future__ import annotations

from concurrent.futures import ProcessPoolExecutor
import functools
import os
import re
import tempfile
from pathlib import Path

import cooler
import pandas as pd
import pyranges1 as pr
from loguru import logger

from capcruncher.types import Assay, BinningMethod
from capcruncher.api.interactions.cooler.merge import merge_coolers

class CoolerBinner:
    def __init__(
        self,
        cooler_group: Path | str | cooler.Cooler,
        binsize: int | None = None,
        method: BinningMethod | str = BinningMethod.MIDPOINT,
        minimum_overlap: float = 0.51,
        n_cis_interaction_correction: bool = True,
        n_rf_per_bin_correction: bool = True,
        scale_factor: int = 1_000_000,
        assay: Assay | str = Assay.CAPTURE,
    ) -> None:
        self.cooler_group = cooler_group
        self.binsize = binsize
        self.method = method
        self.minimum_overlap = minimum_overlap

        if isinstance(cooler_group, str | os.PathLike):
            self.cooler = cooler.Cooler(os.fspath(cooler_group))
        elif isinstance(cooler_group, cooler.Cooler):
            self.cooler = cooler_group
        else:
            raise ValueError(
                "cooler_group must be a path to a cooler file or a cooler object"
            )

        self.n_cis_interactions = self.cooler.info["metadata"]["n_cis_interactions"]
        self.n_cis_interaction_correction = n_cis_interaction_correction
        self.n_restriction_fragment_correction = n_rf_per_bin_correction
        self.scale_factor = scale_factor
        self.assay = assay

    @functools.cached_property
    def genomic_bins(self) -> pr.PyRanges:
        return (
            cooler.binnify(binsize=self.binsize, chromsizes=self.cooler.chromsizes)
            .sort_values(by=["chrom", "start", "end"])
            .assign(
                genomic_bin_id=lambda df: df.reset_index(drop=True)
                .index.to_series()
                .values
            )
            .rename(columns={"chrom": "Chromosome", "start": "Start", "end": "End"})
            .pipe(pr.PyRanges)
        )

    @functools.cached_property
    def fragment_bins(self):
        return (
            self.cooler.bins()[:]
            .rename(
                columns={
                    "chrom": "Chromosome",
                    "start": "Start",
                    "end": "End",
                    "name": "fragment_id",
                }
            )
            .pipe(pr.PyRanges)
        )

    @functools.cached_property
    def fragment_to_genomic_table(self) -> pr.PyRanges:
        """
        Translate genomic bins to fragment bins
        """

        fragment_bins = self.fragment_bins

        if self.method == BinningMethod.MIDPOINT:
            df_fragment_bins = fragment_bins.copy()
            midpoint = (
                df_fragment_bins["Start"].astype(int)
                + (
                    (
                        df_fragment_bins["End"].astype(int)
                        - df_fragment_bins["Start"].astype(int)
                    )
                    // 2
                )
            )
            fragment_bins = pr.PyRanges(
                df_fragment_bins.assign(Start=midpoint, End=midpoint + 1)
            )

        df_fragment_to_bins = self.genomic_bins.join_overlaps(
            fragment_bins,
            strand_behavior="ignore",
            join_type="inner",
            report_overlap_column="Overlap",
        )

        if self.method == BinningMethod.OVERLAP:
            df_fragment_to_bins = df_fragment_to_bins[
                df_fragment_to_bins["Overlap"] >= self.minimum_overlap
            ]

        # Add number of fragments per bin
        df_fragment_to_bins = df_fragment_to_bins.assign(
            n_fragments_per_bin=lambda df: df.groupby("genomic_bin_id")[
                "fragment_id"
            ].transform("nunique"),
        )

        return pr.PyRanges(df_fragment_to_bins)

    @functools.cached_property
    def bins(self) -> pd.DataFrame:
        """Return genomic bins in bedgraph-style column naming."""
        return (
            pd.DataFrame(self.genomic_bins)
            .rename(
                columns={
                    "Chromosome": "chrom",
                    "Start": "start",
                    "End": "end",
                    "genomic_bin_id": "name",
                }
            )[["chrom", "start", "end", "name"]]
            .copy()
        )

    @functools.cached_property
    def bin_conversion_table(self) -> pd.DataFrame:
        """Return fragment-to-genomic-bin mappings using legacy column names."""
        table = pd.DataFrame(self.fragment_to_genomic_table).rename(
            columns={
                "genomic_bin_id": "name_bin",
                "fragment_id": "name_fragment",
                "Overlap": "overlap",
            }
        )
        table["overlap_fraction"] = table["overlap"] / (
            table["End"] - table["Start"]
        )
        return table

    @functools.cached_property
    def fragment_to_genomic_mapping(self) -> dict[int, int]:
        """
        Translate genomic bins to fragment bins
        """
        fragment_to_bins_mapping = (
            self.fragment_to_genomic_table
            .set_index("fragment_id")["genomic_bin_id"]
            .to_dict()
        )
        return fragment_to_bins_mapping

    @functools.cached_property
    def pixels(self) -> pd.DataFrame:
        """
        Translate fragment pixels to genomic pixels
        """

        fragment_to_bins_mapping = self.fragment_to_genomic_mapping

        pixels = self.cooler.pixels()[:].assign(
            genomic_bin1_id=lambda df: df["bin1_id"].map(fragment_to_bins_mapping),
            genomic_bin2_id=lambda df: df["bin2_id"].map(fragment_to_bins_mapping),
        )

        # Sum the counts of pixels that map to the same genomic bins
        pixels = (
            pixels.groupby(["genomic_bin1_id", "genomic_bin2_id"])
            .agg(
                count=("count", "sum"),
            )
            .reset_index()
        )

        # Normalize pixels if specified
        if self.n_restriction_fragment_correction:
            n_fragments_per_bin = (
                self.fragment_to_genomic_table
                .set_index("genomic_bin_id")["n_fragments_per_bin"]
                .to_dict()
            )
            pixels = pixels.assign(
                n_fragments_per_bin1=lambda df: df["genomic_bin1_id"].map(
                    n_fragments_per_bin
                ),
                n_fragments_per_bin2=lambda df: df["genomic_bin2_id"].map(
                    n_fragments_per_bin
                ),
                n_fragments_per_bin_correction=lambda df: (
                    df["n_fragments_per_bin1"] + df["n_fragments_per_bin2"]
                ),
                count_n_rf_norm=lambda df: df["count"]
                / df["n_fragments_per_bin_correction"],
            )

        if self.n_cis_interaction_correction:
            pixels = pixels.assign(
                count_n_cis_norm=lambda df: (df["count"] / self.n_cis_interactions)
                * self.scale_factor,
            )

        if self.n_cis_interaction_correction and self.n_restriction_fragment_correction:
            pixels = pixels.assign(
                count_n_cis_rf_norm=lambda df: (
                    pixels["count_n_rf_norm"] / self.n_cis_interactions
                )
                * self.scale_factor
            )

        return pixels

    @functools.cached_property
    def viewpoint_bins(self) -> list[int]:
        """
        Return list of viewpoint bins
        """

        viewpoint_coords = self.cooler.info["metadata"]["viewpoint_coords"][0]
        match = re.fullmatch(r"([^:]+):(\d+)-(\d+)", viewpoint_coords)
        if match is None:
            raise ValueError(f"Invalid viewpoint coordinates: {viewpoint_coords}")

        chrom, start, end = match.groups()
        pr_viewpoint = pr.PyRanges(
            pd.DataFrame(
                {
                    "Chromosome": [chrom],
                    "Start": [int(start)],
                    "End": [int(end)],
                }
            )
        )

        return pr_viewpoint.join_overlaps(
            self.genomic_bins, strand_behavior="ignore"
        )["genomic_bin_id"].to_list()

    def to_cooler(self, store: Path | str) -> str:
        store = os.fspath(store)
        metadata = {**self.cooler.info["metadata"]}
        metadata["viewpoint_bins"] = [int(x) for x in self.viewpoint_bins]
        metadata["n_interactions_total"] = int(self.cooler.pixels()[:]["count"].sum())
        cooler_fn = f"{store}::/{metadata['viewpoint_name']}/resolutions/{self.binsize}"

        pixels = (
            self.pixels.drop(
                columns=[
                    "bin1_id",
                    "bin2_id",
                    "n_fragments_per_bin1",
                    "n_fragments_per_bin2",
                    "n_fragments_per_bin_correction",
                ],
                errors="ignore",
            )
            .rename(
                columns={"genomic_bin1_id": "bin1_id", "genomic_bin2_id": "bin2_id"}
            )
            .loc[:, lambda df: ["bin1_id", "bin2_id", "count", *df.columns[3:]]]
            .sort_values(by=["bin1_id", "bin2_id"])
        )

        bins = (
            pd.DataFrame(self.genomic_bins).copy()
            .rename(
                columns={"Chromosome": "chrom", "Start": "start", "End": "end"}
            )
            .sort_values("genomic_bin_id")
            .assign(bin_id=lambda df: df["genomic_bin_id"])
            .set_index("genomic_bin_id")
        )

        cooler.create_cooler(
            cooler_fn,
            bins=bins,
            pixels=pixels,
            metadata=metadata,
            mode="w" if not os.path.exists(store) else "a",
            columns=pixels.columns[2:],
            dtypes=dict(zip(pixels.columns[2:], ["float32"] * len(pixels.columns[2:]))),
            ensure_sorted=True,
            ordered=True,
        )

        return cooler_fn


def _bin_cooler(clr_in: str, clr_out: str, binsize: int, **kwargs) -> str:
    clr_binner = CoolerBinner(
        cooler_group=clr_in,
        binsize=binsize,
        **kwargs,
    )
    clr_binner.to_cooler(clr_out)
    return clr_out


def _bin_coolers_local(tasks: list[tuple[str, str, int, dict]]) -> list[str]:
    return [
        _bin_cooler(clr_in, clr_out, binsize, **kwargs)
        for clr_in, clr_out, binsize, kwargs in tasks
    ]


def bins(
    cooler_path: Path | str,
    output: Path | str,
    binsizes: tuple[int, ...] | None = None,
    normalise: bool = True,
    scale_factor: float = 1e6,
    overlap_fraction: float = 1e-9,
    conversion_tables: Path | str | None = None,
    n_cores: int = 1,
    assay: Assay = Assay.CAPTURE,
    **kwargs,
) -> None:
    """
    Convert restriction-fragment cooler groups to constant genomic windows.
    """
    cooler_path = os.fspath(cooler_path)

    clr_groups = cooler.api.list_coolers(cooler_path)

    assert clr_groups, "No cooler groups found in file"
    assert binsizes, "No binsizes provided"

    binning_tasks = []

    for binsize in binsizes:
        for clr_group in clr_groups:
            logger.info(f"Processing {clr_group}")
            clr_in = f"{cooler_path}::{clr_group}"
            clr_out = tempfile.NamedTemporaryFile().name

            default_kwargs = dict(
                method="midpoint",
                minimum_overlap=overlap_fraction,
                n_cis_interaction_correction=normalise,
                n_rf_per_bin_correction=normalise,
                scale_factor=scale_factor,
                assay=assay,
            )

            binning_tasks.append((clr_in, clr_out, binsize, default_kwargs | kwargs))

    if n_cores > 1 and len(binning_tasks) > 1:
        try:
            with ProcessPoolExecutor(max_workers=n_cores) as executor:
                futures = [
                    executor.submit(_bin_cooler, clr_in, clr_out, binsize, **kwargs)
                    for clr_in, clr_out, binsize, kwargs in binning_tasks
                ]
                clr_tempfiles = [future.result() for future in futures]
        except OSError as exc:
            logger.warning(
                f"Process executor unavailable ({exc}); falling back to local binning."
            )
            clr_tempfiles = _bin_coolers_local(binning_tasks)
    else:
        clr_tempfiles = _bin_coolers_local(binning_tasks)

    merge_coolers([Path(clr_tempfile) for clr_tempfile in clr_tempfiles], output)
