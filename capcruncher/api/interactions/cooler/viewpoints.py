from __future__ import annotations

from pathlib import Path
from typing import Self

import pyranges1 as pr

from capcruncher.types import Assay, VALID_ASSAYS, validate_choice

class Viewpoint:
    def __init__(self, coordinates: pr.PyRanges, assay: Assay | str) -> None:
        self.coordinates = coordinates
        self.assay = validate_choice(assay, VALID_ASSAYS, "assay")

    @classmethod
    def from_bed(
        cls, bed: Path | str, viewpoint: str, assay: Assay | str
    ) -> Self:
        """
        Creates a viewpoint object from a bed file.

        Args:
            bed (str): Path to bed file containing viewpoint coordinates.
            viewpoint (str): Name of viewpoint to extract from bed file.

        Raises:
            IndexError: Oligo name cannot be found within viewpoints.

        Returns:
            Viewpoint: Viewpoint object.
        """
        df_viewpoints = pr.read_bed(Path(bed))
        viewpoint_names = df_viewpoints["Name"].astype(str)
        df_viewpoints = df_viewpoints.loc[
            (viewpoint_names == viewpoint)
            | viewpoint_names.str.endswith(f"_{viewpoint}")
        ]

        if df_viewpoints.empty:
            raise IndexError(
                f"Oligo name cannot be found within viewpoints: {viewpoint}"
            )

        return cls(pr.PyRanges(df_viewpoints), assay=assay)

    def bins(self, bins: pr.PyRanges):
        """
        Returns the bins that overlap with the viewpoint.

        Args:
            bins (pr.PyRanges): PyRanges object containing all bins.

        Returns:
            pr.PyRanges: PyRanges object containing all bins that overlap with the viewpoint.
        """
        return bins.join_overlaps(self.coordinates, strand_behavior="ignore")

    def bin_names(self, bins: pr.PyRanges) -> list[int]:
        return self.bins(bins)["Name"].astype(int).to_list()

    def bins_cis(self, bins: pr.PyRanges) -> list[int]:
        """
        Returns the bins that are on the same chromosome(s) as the viewpoint.

        Args:
            bins (pr.PyRanges): PyRanges object containing all bins.

        Returns:
            List[int]: List of bin names.
        """

        # Get the chromosomes of the viewpoint
        viewpoint_chromosomes = self.chromosomes

        # Get the bins that are on the same chromosome(s) as the viewpoint
        df_cis_bins = bins.loc[
            lambda df: df["Chromosome"].isin(viewpoint_chromosomes)
        ]

        # If capture or tri, remove viewpoint bins from cis bins
        if self.assay in {Assay.CAPTURE, Assay.TRI}:
            df_cis_bins = df_cis_bins.loc[
                lambda df: ~df["Name"].isin(self.bin_names(bins))
            ]

        return df_cis_bins["Name"].astype(int).to_list()

    @property
    def chromosomes(self) -> list[str]:
        return self.coordinates["Chromosome"].unique().tolist()

    @property
    def coords(self) -> list[str]:
        """
        Returns the genomic coordinates of the viewpoint.

        Returns:
            List[str]: List of genomic coordinates.
        """
        _coords = []
        for row in self.coordinates.itertuples():
            _coords.append(f"{row.Chromosome}:{row.Start}-{row.End}")

        return _coords

