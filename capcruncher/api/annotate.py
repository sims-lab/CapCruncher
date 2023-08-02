import warnings
from typing import Union

import pandas as pd
import pybedtools
import pyranges as pr
import ray
from loguru import logger

from capcruncher.utils import convert_bed_to_pr

warnings.simplefilter("ignore", category=RuntimeWarning)


@ray.remote
class BedFileIntersection:
    def __init__(
        self,
        bed_a: Union[str, pybedtools.BedTool, pr.PyRanges],
        bed_b: Union[str, pybedtools.BedTool, pr.PyRanges],
        name: str = "b",
        action: str = "get",
        fraction: float = 1e-9,
    ):

        self.a = bed_a
        self.b = bed_b
        self.name = name
        self.action = action
        self.fraction = fraction

        self.pr_a = convert_bed_to_pr(self.a, ignore_ray_objrefs=True)

        import logging

        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

    def _get_intersection(self, pr_b: pr.PyRanges):

        intersection = (
            self.pr_a.join(
                pr_b,
                report_overlap=True,
            )
            .assign("frac", lambda df: df.eval("Overlap / (End - Start)"))
            .subset(lambda df: df["frac"] >= self.fraction)
            .as_df()
        )

        dtype = pd.CategoricalDtype(pr_b.df["Name"].unique())

        intersection_data = (
            intersection.set_index("Name")["Name_b"].astype(dtype).rename(self.name)
        )

        return intersection_data

    def _count_intersection(self, pr_b: pr.PyRanges):

        intersection_data = (
            self.pr_a.coverage(pr_b)
            .df.query(f"NumberOverlaps > 0 and FractionOverlaps >= {self.fraction}")
            .set_index("Name")["NumberOverlaps"]
            .rename(self.name)
        )

        return intersection_data

    def intersection(self):

        try:

            pr_b = convert_bed_to_pr(self.b)

            if self.action == "get":
                _intersection = self._get_intersection(pr_b)
            elif self.action == "count":
                _intersection = self._count_intersection(pr_b)

        except (OSError, IndexError, FileNotFoundError, StopIteration, AssertionError):

            self.logger.warning(
                f"Could not intersect {self.b} using {self.action} method."
            )
            _intersection = pd.Series(
                data=pd.NA,
                index=self.pr_a.df["Name"],
                name=self.name,
                dtype=object,
            )

        return _intersection

    def __repr__(self):
        return f"{self.name} intersection"
