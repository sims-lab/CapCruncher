import itertools
from typing import Union
import warnings

warnings.simplefilter("ignore")

import pandas as pd
from joblib import Parallel, delayed
from pybedtools import BedTool

from ccanalyser.utils import (
    bed_has_name,
    convert_bed_to_dataframe,
    convert_to_bedtool,
    is_valid_bed,
    split_intervals_on_chrom,
)
from natsort import natsorted
import traceback


class BedIntersection:
    def __init__(
        self,
        bed1: Union[str, BedTool, pd.DataFrame],
        bed2: Union[str, BedTool, pd.DataFrame],
        intersection_name: str = "count",
        intersection_method: str = "count",
        intersection_min_frac: float = 1e-9,
        n_cores: int = 1,
        invalid_bed_action="error",
    ):

        self.bed1 = bed1
        self.bed2 = bed2
        self.methods = {
            "get": self._intersections_get,
            "count": self._intersections_count,
        }
        self.methods_na = {"get": ".", "count": 0}
        self.intersection_name = intersection_name
        self.intersection_method = self.methods.get(intersection_method)
        self.intersection_na = self.methods_na.get(intersection_method)
        self.min_frac = intersection_min_frac

        self.invalid_bed_action = invalid_bed_action
        self.bed1_valid = is_valid_bed(bed1)
        self.bed2_valid = is_valid_bed(bed2)

        self.n_cores = n_cores

    def _intersections_count(self, a, b):
        return a.intersect(
            b, loj=True, c=True, f=self.min_frac, sorted=True
        ).to_dataframe()

    def _intersections_get(self, a, b):
        return a.intersect(b, loj=True, f=self.min_frac, sorted=True).to_dataframe()

    def _extract_series(self, intersection):
        return intersection.to_dataframe().iloc[:, -1]

    def _intersect_by_chrom(
        self, a: Union[BedTool, pd.DataFrame], b: Union[BedTool, pd.DataFrame]
    ):

        a_by_chrom = split_intervals_on_chrom(a)
        b_by_chrom = split_intervals_on_chrom(b)

        a_chroms, b_chroms = set(a_by_chrom), set(b_by_chrom)

        intersection_required = []
        not_intersected = []
        for chrom in a_chroms:

            a_chrom = a_by_chrom[chrom]

            if chrom in b_chroms:
                b_chrom = b_by_chrom[chrom]

                a_chrom_bed = convert_to_bedtool(a_chrom)
                b_chrom_bed = convert_to_bedtool(b_chrom)

                intersection_required.append(
                    delayed(self.intersection_method)(a_chrom_bed, b_chrom_bed)
                )
            else:
                not_intersected.append(a_chrom)

        intersections = Parallel(n_jobs=self.n_cores)(intersection_required)

        return pd.concat([*intersections, *not_intersected], ignore_index=True).fillna(
            self.intersection_na
        )

    @property
    def intersection(self) -> pd.Series:

        if all([self.bed1_valid, self.bed2_valid]):
            a = convert_to_bedtool(self.bed1)
            b = convert_to_bedtool(self.bed2)
            _intersection = self._intersect_by_chrom(a, b)
            ser = _intersection.set_index("name").iloc[:, -1]
            ser.name = self.intersection_name

        elif self.invalid_bed_action == "ignore":
            ser = (
                convert_bed_to_dataframe(self.bed1)
                .assign(**{self.intersection_name: self.intersection_na})
                .set_index("name")
                .loc[:, self.intersection_name]
            )
        else:
            raise ValueError(
                f"Slices valid: {self.bed1_valid}\n {self.intersection_name} .bed file valid: {self.bed2_valid}"
            )
            traceback.format_exc()

        return ser