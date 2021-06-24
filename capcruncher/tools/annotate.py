from typing import Union
import warnings

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
        intersection_split_chrom: bool = True, 
        n_cores: int = 1,
        invalid_bed_action="error",
        
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
        self.intersection_split_chrom = intersection_split_chrom
        self.n_cores = n_cores
        self.invalid_bed_action = invalid_bed_action

    def _intersections_count(self, a, b):
        return a.intersect(
            b, loj=True, c=True, f=self.min_frac, sorted=True).to_dataframe()

    def _intersections_get(self, a, b):
        return a.intersect(b, loj=True, f=self.min_frac, sorted=True).to_dataframe()

    def _extract_series(self, intersection):
        return intersection.to_dataframe().iloc[:, -1]

    def _intersect_by_chrom(
        self, a: Union[BedTool, pd.DataFrame], b: Union[BedTool, pd.DataFrame]
    ):
        
        from joblib import Parallel, delayed

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
                    delayed(self._intersection_method)(a_chrom_bed, b_chrom_bed)
                )
            else:
                not_intersected.append(a_chrom)

        intersections = Parallel(n_jobs=self.n_cores)(intersection_required)

        return pd.concat([*intersections, *not_intersected], ignore_index=True).fillna(
            self._intersection_na
        )
    

    def _format_invalid_intersection(self, bed):
        return (
                convert_bed_to_dataframe(bed)
                .assign(**{self.intersection_name: self._intersection_na})
                .set_index("name")
                .loc[:, self.intersection_name])


    @property
    def intersection(self) -> pd.Series:
        '''Intersects the two bed files and returns a pd.Series.'''

        if all([self.bed1_valid, self.bed2_valid]):

            a = convert_to_bedtool(self.bed1)
            b = convert_to_bedtool(self.bed2)
            
            if self.intersection_split_chrom:
                _intersection = self._intersect_by_chrom(a, b)
            else:
                _intersection = self._intersection_method(a, b)
            
            ser = _intersection.set_index("name").iloc[:, -1]
            ser.name = self.intersection_name


        elif self.invalid_bed_action == "ignore":
            ser = self._format_invalid_intersection(self.bed1)
        
        else:
            raise ValueError(
                f"Slices valid: {self.bed1_valid}\n {self.intersection_name} .bed file valid: {self.bed2_valid}"
            )     

        return ser