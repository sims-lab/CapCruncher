from logging import disable
import logging
from typing import Union, Literal
import warnings
import numpy as np

warnings.simplefilter("ignore", category=RuntimeWarning)

import pandas as pd
from pybedtools import BedTool

from capcruncher.utils import (
    convert_bed_to_dataframe,
    convert_to_bedtool,
    is_valid_bed,
    split_intervals_on_chrom,
)


class BedIntersection:

    """
    Performs intersection between two named bed files.

    Wrapper around bedtools intersect designed to intersect in parallel
    (by splitting file based on chromosome) and handle malformed bed files.

    Attributes:
     bed1 (Union[str, BedTool, pd.DataFrame]): Bed file to intersect. Must be named.
     bed2 (Union[str, BedTool, pd.DataFrame]): Bed file to intersect.
     intersection_name (str): Name for intersection.
     min_frac (float): Minimum fraction required for intersection
     n_cores (int): Number of cores for parallel intersection.
     invalid_bed_action: Method to deal with missing/malformed bed files ("ignore"|"error")

    """

    def __init__(
        self,
        bed1: Union[str, BedTool, pd.DataFrame],
        bed2: Union[str, BedTool, pd.DataFrame],
        intersection_name: str = "count",
        intersection_method: str = "count",
        intersection_min_frac: float = 1e-9,
        invalid_bed_action: Literal["ignore", "error"] = "error",
        dtype: str = "str",
    ):
        """

        Args:
         bed1 (Union[str, BedTool, pd.DataFrame]): Bed file 1 for intersection
         bed2 (Union[str, BedTool, pd.DataFrame]): Bed file 2 for intersection
         intersection_name (str, optional): Name for intersection. Defaults to "count".
         intersection_method (str, optional): Method for intersection. Choose from (get|count). Defaults to "count".
         intersection_min_frac (float, optional): Minimum overlap fraction required. Defaults to 1e-9.
         n_cores (int, optional): Number of cores to use for intersection. Defaults to 1.
         invalid_bed_action (str, optional): If set to 'ignore' will return default N/A value or raise Exception. Defaults to "error".
        """

        # Format input bed files
        self.bed1 = bed1
        self.bed2 = bed2

        self.bed1_valid = is_valid_bed(bed1)
        self.bed2_valid = is_valid_bed(bed2)

        # Default methods and NA values
        self._methods = {
            "get": self._intersections_get,
            "count": self._intersections_count,
        }

        self._methods_na = {"get": ".", "count": 0}

        # Get intersection properties
        self.intersection_name = intersection_name
        self._intersection_method = self._methods.get(intersection_method)
        self._intersection_na = self._methods_na[intersection_method]
        self.min_frac = intersection_min_frac

        # Other options
        self.invalid_bed_action = invalid_bed_action
        self.dtype = dtype

    def _extract_intersection_result(
        self, bt: BedTool, intersection_id: str, dtype: str = None
    ):
        columns = pd.read_csv(bt.fn, nrows=5, sep="\t", header=None).columns
        columns_needed = [columns[3], columns[-1]]

        return pd.read_csv(
            bt.fn,
            sep="\t",
            header=None,
            usecols=columns_needed,
            names=["parent_name", intersection_id],
            index_col="parent_name",
            squeeze=True,
            dtype={"parent_name": str, intersection_id: dtype},
            na_values=".",
        )

    def _intersections_count(self, a, b):
        return self._extract_intersection_result(
            a.intersect(b, loj=True, c=True, f=self.min_frac, sorted=True),
            intersection_id=self.intersection_name,
            dtype=self.dtype,
        )

    def _intersections_get(self, a, b):

        if self.dtype == "category":
            b_names = b.to_dataframe().iloc[:, -1].unique()
            dtype = pd.CategoricalDtype(b_names)
        else:
            dtype = self.dtype

        return self._extract_intersection_result(
            a.intersect(b, loj=True, f=self.min_frac, sorted=True),
            intersection_id=self.intersection_name,
            dtype=dtype,
        )

    def _format_invalid_intersection(self, bed):
        return (
            convert_bed_to_dataframe(bed)
            .assign(**{self.intersection_name: np.nan})
            .set_index("name")
            .loc[:, self.intersection_name]
            .astype(float)
        )

    def get_intersection(self, sort_index: bool = True) -> pd.Series:
        """Intersects the two bed files and returns a pd.Series."""

        logging.info(f"Performing {self.intersection_name} intersection.")

        if all([self.bed1_valid, self.bed2_valid]):

            a = convert_to_bedtool(self.bed1)
            b = convert_to_bedtool(self.bed2)
            intersection = self._intersection_method(a, b)

        elif self.invalid_bed_action == "ignore":
            intersection = self._format_invalid_intersection(self.bed1)

        else:
            raise ValueError(
                f"Slices valid: {self.bed1_valid}\n {self.intersection_name} .bed file valid: {self.bed2_valid}"
            )

        return intersection.sort_index() if sort_index else intersection
