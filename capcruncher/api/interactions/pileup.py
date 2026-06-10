from __future__ import annotations

import os
import subprocess
import tempfile
from pathlib import Path

import cooler
from loguru import logger
from pydantic import BaseModel, Field, PositiveFloat, field_validator, model_validator

from capcruncher.api.interactions.bedgraph import CoolerBedGraph
from capcruncher.types import Normalisation, PileupFormat


class PileupOptions(BaseModel):
    """Validated options for extracting bedgraph or bigWig pileups."""

    uri: Path | str
    viewpoint_names: list[str] | None = None
    output_prefix: Path | str = ""
    format: PileupFormat = PileupFormat.BEDGRAPH
    normalisation: Normalisation = Normalisation.RAW
    normalisation_regions: Path | str | None = None
    binsize: int = Field(default=0, ge=0)
    gzip: bool = True
    scale_factor: PositiveFloat = 1e6
    sparse: bool = True

    @field_validator("normalisation_regions", mode="before")
    @classmethod
    def empty_region_to_none(cls, value: Path | str | None) -> Path | str | None:
        return None if value == "" else value

    @model_validator(mode="after")
    def validate_normalisation_regions(self) -> PileupOptions:
        if (
            self.normalisation == Normalisation.REGION
            and self.normalisation_regions is None
        ):
            raise ValueError(
                "normalisation_regions is required when normalisation is 'region'."
            )
        if (
            self.normalisation != Normalisation.REGION
            and self.normalisation_regions is not None
        ):
            raise ValueError(
                "normalisation_regions can only be used when normalisation is 'region'."
            )
        return self


def pileup(
    uri: Path | str,
    viewpoint_names: list[str] | None = None,
    output_prefix: Path | str = "",
    format: PileupFormat = PileupFormat.BEDGRAPH,
    normalisation: Normalisation = Normalisation.RAW,
    normalisation_regions: Path | str | None = None,
    binsize: int = 0,
    gzip: bool = True,
    scale_factor: float = 1e6,
    sparse: bool = True,
) -> None:
    """Extract reporters from a capture experiment as bedgraph or bigWig files.

    Identifies reporters for one viewpoint, if supplied, or all capture probes present
    in a CapCruncher HDF5 file.
    """
    options = PileupOptions(
        uri=uri,
        viewpoint_names=viewpoint_names,
        output_prefix=output_prefix,
        format=format,
        normalisation=normalisation,
        normalisation_regions=normalisation_regions,
        binsize=binsize,
        gzip=gzip,
        scale_factor=scale_factor,
        sparse=sparse,
    )

    uri = os.fspath(options.uri)
    output_prefix = os.fspath(options.output_prefix)
    normalisation_regions = (
        os.fspath(options.normalisation_regions)
        if options.normalisation_regions is not None
        else None
    )
    viewpoint_names = options.viewpoint_names or [
        v.strip("/") for v in cooler.fileops.list_coolers(uri) if "resolutions" not in v
    ]

    logger.info(f"Performing pileup for {viewpoint_names}")

    bin_bedgraph = options.binsize > 0

    for viewpoint_name in viewpoint_names:
        cooler_group = f"{uri}::{viewpoint_name}"

        if bin_bedgraph:
            cooler_group = f"{cooler_group}/resolutions/{options.binsize}"

        try:
            cooler.fileops.is_cooler(cooler_group)
        except Exception as exc:
            logger.warning(
                f"Viewpoint {viewpoint_name} not found in cooler "
                f"(cooler may be empty): {exc}. Writing empty output."
            )
            if options.format == PileupFormat.BEDGRAPH:
                open(
                    f"{output_prefix}_{viewpoint_name}.bedgraph{'.gz' if options.gzip else ''}",
                    "w",
                ).close()
            continue

        bedgraph = CoolerBedGraph(uri=cooler_group, sparse=sparse).extract_bedgraph(
            normalisation=options.normalisation,
            region=normalisation_regions,
            scale_factor=options.scale_factor,
        )

        logger.info(f"Generated bedgraph for {viewpoint_name}")

        if options.format == PileupFormat.BEDGRAPH:
            bedgraph.to_csv(
                f"{output_prefix}_{viewpoint_name}.bedgraph{'.gz' if options.gzip else ''}",
                sep="\t",
                header=False,
                index=False,
            )
        elif options.format == PileupFormat.BIGWIG:
            clr = cooler.Cooler(cooler_group)

            with tempfile.NamedTemporaryFile() as chromsizes_tmp:
                with tempfile.NamedTemporaryFile() as bedgraph_tmp:
                    clr.chromsizes.to_csv(chromsizes_tmp, sep="\t", header=False)
                    bedgraph.to_csv(bedgraph_tmp, sep="\t", index=False, header=False)

                    subprocess.run(
                        [
                            "bedGraphToBigWig",
                            bedgraph_tmp.name,
                            chromsizes_tmp.name,
                            f"{output_prefix}_{viewpoint_name}.bigWig",
                        ],
                        check=True,
                    )
