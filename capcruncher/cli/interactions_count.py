import os
from loguru import logger
import tempfile
import glob
from typing import Literal


def count(
    reporters: os.PathLike,
    output: os.PathLike = "CC_cooler.hdf5",
    remove_exclusions: bool = False,
    remove_viewpoint: bool = False,
    subsample: float = 0,
    fragment_map: os.PathLike = None,
    viewpoint_path: os.PathLike = None,
    n_cores: int = 1,
    assay: Literal["capture", "tri", "tiled"] = "capture",
    **kwargs,
) -> os.PathLike:
    """
    Counts interactions between the viewpoint and the rest of the genome.

    Args:
        reporters: Path to reporters file.
        output: Output file name.
        remove_exclusions: Remove excluded regions.
        remove_viewpoint: Remove capture regions.
        subsample: Subsample reads.
        fragment_map: Path to fragment map.
        viewpoint_path: Path to viewpoint file.
        n_cores: Number of cores.
        assay: Assay type.
        **kwargs: Additional arguments.
    Returns:
        Path to the generated cooler file.

    """
    from capcruncher_tools.api import count_interactions

    clr = count_interactions(
        reporters=reporters,
        output=output,
        remove_exclusions=remove_exclusions,
        remove_viewpoint=remove_viewpoint,
        subsample=subsample,
        fragment_map=fragment_map,
        viewpoint_path=viewpoint_path,
        n_cores=n_cores,
        assay=assay,
        **kwargs,
    )

    return clr
