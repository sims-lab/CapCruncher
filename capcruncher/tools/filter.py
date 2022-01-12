import logging
import pandas as pd
import os
import numpy as np
import itertools


class SliceFilter:

    """
    Perform slice filtering (inplace) and reporter identification.

    The SliceFilter classes e.g. CCSliceFilter, TriCSliceFilter, TiledCSliceFilter
    perform all of the filtering (inplace) and reporter identification whilst also
    providing statistics of the numbers of slices/reads removed at each stage.

    Attributes:
     slices (pd.DataFrame): Annotated slices dataframe.
     fragments (pd.DataFrame): Slices dataframe aggregated by parental read.
     reporters (pd.DataFrame): Slices identified as reporters.
     filter_stages (dict): Dictionary containg stages and a list of class methods (str) required to get to this stage.
     slice_stats (pd.DataFrame): Provides slice level statistics.
     read_stats (pd.DataFrame): Provides statistics of slice filtering at the parental read level.
     filter_stats (pd.DataFrame): Provides statistics of read filtering.

    """

    def __init__(
        self,
        slices: pd.DataFrame,
        filter_stages: dict = None,
        sample_name: str = "",
        read_type: str = "",
    ):
        """
        Base for all slice filter objects.

        Slices DataFrame must have the following columns:

         - slice_name: Unique aligned read identifier (e.g. XZKG:889:11|flashed|1)
         - parent_read: Identifier shared by slices from same fragment (e.g.XZKG:889:11)
         - pe: Read combined by FLASh or not (i.e. "flashed" or "pe")
         - mapped: Alignment is mapped (e.g. 0/1)
         - multimapped: Alignment is mapped (e.g. 0/1)
         - slice: Slice number (e.g. 0)
         - chrom: Chromosome e.g. chr1
         - start: Start coord
         - end: End coord
         - capture: Capture site intersecting slice (e.g. Slc25A37)
         - capture_count: Number of capture probes overlapping slice (e.g. 1)
         - exclusion: Read present in excluded region (e.g. Slc25A37)
         - exclusion_count: Number of excluded regions overlapping slice (e.g. 1)
         - blacklist: Read present in excluded region (e.g. 0)
         - coordinates: Genome coordinates (e.g. chr1:1000-2000)

        Filtering to be performed can be left as the default (all start with 'remove')
        or a custom filtering order can be supplied with a yaml file. This must have the format:

         FILTER_STAGE_NAME:
             - FILTER 1
             - FILTER 2
         FILTER_STAGE_NAME2:
             - FILTER 3
             - FILTER 1


        *All* filters present in the file must be defined within the SliceFilter class.


        Args:
         slices (pd.DataFrame): DatFrame containing annotated slices
         filter_stages (dict, optional): Dictionary defining order of slice filtering. Defaults to None.
         sample_name (str, optional): Name of sample being processed e.g. DOX-treated_1. Defaults to "".
         read_type (str, optional): Combined (flashed) or not-combined (pe). Defaults to "".

        Raises:
         ValueError: Filter stages must be provided. This is done automatically by all subclasses
         AttributeError: All filters must be defined in the SliceFilter.
        """

        self._has_required_columns = self._required_columns_present(slices)
        self.slices = slices.sort_values(["parent_read", "slice"])

        if filter_stages:
            self.filter_stages = self._extract_filter_stages(filter_stages)
        else:
            raise ValueError("Filter stages not provided")

        self._filter_stats = pd.DataFrame()
        self.sample_name = sample_name
        self.read_type = read_type

    def _required_columns_present(self, df) -> bool:

        columns_required = [
            "parent_id",
            "slice_name",
            "parent_read",
            "pe",
            "mapped",
            "multimapped",
            "slice",
            "chrom",
            "start",
            "end",
            "capture",
            "capture_count",
            "exclusion",
            "blacklist",
            "coordinates",
        ]

        for col in columns_required:
            if not col in df.columns:
                raise KeyError(f'Required column "{col}" not in slices dataframe')

        return True

    def _extract_filter_stages(self, filter_stages) -> dict:
        """
        Extracts filter stages from a supplied dictionary or yaml file

        Checks that the filters provided are within the dictionary supplied.
        """

        if isinstance(filter_stages, dict):
            filters = filter_stages

        elif os.path.exists(filter_stages) and (
            ".yaml" in filter_stages or ".yml" in filter_stages
        ):

            import yaml

            with open(filter_stages, "r") as f:
                filters = yaml.safe_load(f)

        else:
            raise ValueError(
                "Provide either a path to a .yaml file or a python dictionary"
            )

        all_filters = itertools.chain.from_iterable(filters.values())

        for filt in all_filters:
            if not filt in self.filters:
                raise AttributeError(
                    f"Required filter: {filt} not present. Check for correct spelling and format."
                )

        return filters

    @property
    def filters(self) -> list:
        """A list of the callable filters present within the slice filterer instance.

        Returns:
            list: All filters present in the class.
        """
        filters = [attr for attr in dir(self) if "remove_" in attr]

        # There is at least one filter not indicated by remove
        # Need to append to the filter list.
        filters.append("get_unfiltered_slices")

        return filters

    @property
    def slice_stats(self) -> pd.DataFrame:
        """
        Statistics at the slice level.

        Returns:
         pd.DataFrame: Statistics per slice.
        """
        raise NotImplementedError("Override this method")

    @property
    def filter_stats(self) -> pd.DataFrame:
        """
        Statistics for each filter stage.

        Returns:
         pd.DataFrame: Statistics of the number of slices removed at each stage.
        """
        return (
            self._filter_stats.transpose()
            .reset_index()
            .rename(columns={"index": "stage"})
            .assign(sample=self.sample_name, read_type=self.read_type)
        )

    @property
    def read_stats(self) -> pd.DataFrame:
        """
        Gets statistics at a read level.

        Aggregates slices by parental read id and calculates stats.

        Returns:
         pd.DataFrame: Statistics of the slices/fragments removed aggregated by read id.
        """
        return self.filter_stats.rename(
            columns={
                "stage": "stat_type",
                "unique_fragments": "stat",
            }
        )[["stat_type", "stat"]].assign(
            stage="ccanalysis",
            read_type=self.read_type,
            sample=self.sample_name,
            read_number=0,
        )

    @property
    def fragments(self) -> pd.DataFrame:
        """
        Summarises slices at the fragment level.

        Uses pandas groupby to aggregate slices by their parental read name
        (shared by all slices from the same fragment). Also determines the
        number of reporter slices for each fragment.

        Returns:
         pd.DataFrame: Slices aggregated by parental read name.

        """
        raise NotImplementedError("Override this property")

    @property
    def captures(self) -> pd.DataFrame:
        raise NotImplementedError("Override this property")

    @property
    def reporters(self) -> pd.DataFrame:
        """
        Extracts reporter slices from slices dataframe i.e. non-capture slices

        Returns:
         pd.DataFrame: All non-capture slices

        """
        raise NotImplementedError("Override this property")

    def filter_slices(self, output_slices=False, output_location="."):
        """
        Performs slice filtering.

        Filters are applied to the slices dataframe in the order specified by
        filter_stages. Filtering stats aggregated at the slice and fragment level
        are also printed.

        Args:
         output_slices (bool, optional): Determines if slices are to be output to a specified location after each filtering step.
                                         Useful for debugging. Defaults to False.
         output_location (str, optional): Location to output slices at each stage. Defaults to ".".
        """

        for stage, filters in self.filter_stages.items():
            for filt in filters:
                # Call all of the filters in the filter_stages dict in order
                logging.info(f"Filtering slices: {filt}")
                getattr(self, filt)()  # Gets and calls the selected method
                logging.info(f"Completed: {filt}")
                logging.info(f"Number of slices: {self.slices.shape[0]}")
                logging.info(f'Number of reads: {self.slices["parent_read"].nunique()}')

                if output_slices == "filter":
                    self.slices.to_csv(os.path.join(output_location, f"{filt}.tsv.gz"))

            if output_slices == "stage":
                self.slices.to_csv(os.path.join(output_location, f"{stage}.tsv.gz"))

            self._filter_stats[stage] = self.slice_stats

    def get_unfiltered_slices(self):
        """
        Does not modify slices.
        """
        self.slices = self.slices

    def remove_unmapped_slices(self):
        """
        Removes slices marked as unmapped (Uncommon)
        """
        self.slices = self.slices.query("mapped == 1")

    def remove_orphan_slices(self):
        """Remove fragments with only one aligned slice (Common)"""

        # fragments = self.fragments
        # fragments_multislice = fragments.query("unique_slices > 1")
        # self.slices = self.slices[
        #     self.slices["parent_read"].isin(fragments_multislice["parent_read"])
        # ]

        not_orphan = self.slices["parent_id"].duplicated(keep=False)
        self.slices = self.slices.loc[not_orphan]

    def remove_duplicate_re_frags(self):
        """
        Prevent the same restriction fragment being counted more than once (Uncommon).

        Example:

         --RE_FRAG1--\----Capture----\---RE_FRAG1----

        """
        self.slices = self.slices.drop_duplicates(
            subset=["parent_read", "restriction_fragment"]
        )

    def remove_slices_without_re_frag_assigned(self):
        """Removes slices if restriction_fragment column is N/A"""
        self.slices = self.slices.query('restriction_fragment != "."')

    def remove_duplicate_slices(self):
        """
        Remove all slices if the slice coordinates and slice order are shared.

        This method is designed to remove a fragment if it is a PCR duplicate
        (Common).

        Example:

         | Frag 1:  chr1:1000-1250 chr1:1500-1750
         | Frag 2:  chr1:1000-1250 chr1:1500-1750
         | Frag 3:  chr1:1050-1275 chr1:1600-1755
         | Frag 4:  chr1:1500-1750 chr1:1000-1250

         Frag 2 removed. Frag 1,3,4 retained


        """

        frags_deduplicated = (
            self.slices.groupby("parent_id")
            .agg(coords=("coordinates", "|".join))
            .reset_index()
            .drop_duplicates(subset="coords", keep="first")
        )

        self.slices = self.slices.loc[
            self.slices["parent_id"].isin(frags_deduplicated["parent_id"])
        ]

    def remove_duplicate_slices_pe(self):
        """
        Removes PCR duplicates from non-flashed (PE) fragments (Common).

        Sequence quality is often lower at the 3' end of reads leading to variance
        in mapping coordinates.  PCR duplicates are removed by checking that the
        fragment start and end are not duplicated in the dataframe.

        """
        if (
            self.slices["pe"].iloc[:100].str.contains("pe").sum() > 1
        ):  # at least one un-flashed

            fragments_partial = (
                self.slices.groupby("parent_id")
                .agg(coords=("coordinates", "|".join))
                .reset_index()
            )

            fragments_partial = fragments_partial.assign(
                read_start=lambda df: df["coords"]
                .str.split("|")
                .str[0]
                .str.split(r":|-")
                .str[1],
                read_end=lambda df: df["coords"]
                .str.split("|")
                .str[-1]
                .str.split(r":|-")
                .str[-1],
            )

            fragments_deduplicated = fragments_partial.drop_duplicates(
                subset=["read_start", "read_end"]
            )

            self.slices = (
                self.slices.set_index("parent_id")
                .loc[fragments_deduplicated["parent_id"]]
                .reset_index()
            )

    def remove_excluded_slices(self):
        """Removes any slices in the exclusion region (default 1kb) (V. Common)"""

        slices_with_viewpoint = self.slices_with_viewpoint
        slices_passed = slices_with_viewpoint.loc[lambda df: (df["exclusion_count"] < 1) | (df["exclusion"] != df["viewpoint"])] 
        
        self.slices = self.slices.loc[
            lambda df: df["parent_id"].isin(slices_passed["parent_id"])
        ]

    def remove_blacklisted_slices(self):
        """Removes slices marked as being within blacklisted regions"""
        self.slices = self.slices.loc[lambda df: (df["blacklist"] == 0) | (df["blacklist"].isna())]

    @property
    def slices_with_viewpoint(self):

        slices = self.slices.set_index("parent_id")
        captures = self.captures.set_index("parent_id")
        return (
            slices.join(captures["capture"], lsuffix="_slices", rsuffix="_capture")
            .rename(
                columns={"capture_slices": "capture", "capture_capture": "viewpoint"}
            )
            .reset_index()
        )


class CCSliceFilter(SliceFilter):
    """
    Perform Capture-C slice filtering (inplace) and reporter identification.

    SliceFilter tuned specifically for Capture-C data. This class has addtional methods
    to remove common artifacts in Capture-C data i.e. multi-capture fragments,
    non-reporter fragments, multi-capture reporters. The default filter order is as follows:

     - remove_unmapped_slices
     - remove_orphan_slices
     - remove_multi_capture_fragments
     - remove_excluded_slices
     - remove_blacklisted_slices
     - remove_non_reporter_fragments
     - remove_viewpoint_adjacent_restriction_fragments
     - remove_slices_without_re_frag_assigned
     - remove_duplicate_re_frags
     - remove_duplicate_slices
     - remove_duplicate_slices_pe
     - remove_non_reporter_fragments

    See the individual methods for further details.

    Attributes:
     slices (pd.DataFrame): Annotated slices dataframe.
     fragments (pd.DataFrame): Slices dataframe aggregated by parental read.
     reporters (pd.DataFrame): Slices identified as reporters.
     filter_stages (dict): Dictionary containg stages and a list of class methods (str) required to get to this stage.
     slice_stats (pd.DataFrame): Provides slice level statistics.
     read_stats (pd.DataFrame): Provides statistics of slice filtering at the parental read level.
     filter_stats (pd.DataFrame): Provides statistics of read filtering.

    """

    def __init__(self, slices, filter_stages=None, **sample_kwargs):
        if not filter_stages:
            filter_stages = {
                "pre-filtering": [
                    "get_unfiltered_slices",
                ],
                "mapped": [
                    "remove_unmapped_slices",
                ],
                "contains_single_capture": [
                    "remove_orphan_slices",
                    "remove_multi_capture_fragments",
                ],
                "contains_capture_and_reporter": [
                    "remove_excluded_slices",
                    "remove_blacklisted_slices",
                    "remove_non_reporter_fragments",
                    "remove_viewpoint_adjacent_restriction_fragments",
                ],
                "duplicate_filtered": [
                    "remove_slices_without_re_frag_assigned",
                    "remove_duplicate_re_frags",
                    "remove_duplicate_slices",
                    "remove_duplicate_slices_pe",
                    "remove_non_reporter_fragments",
                ],
            }

        super(CCSliceFilter, self).__init__(slices, filter_stages, **sample_kwargs)

    @property
    def fragments(self) -> pd.DataFrame:
        """
        Summarises slices at the fragment level.

        Uses pandas groupby to aggregate slices by their parental read name
        (shared by all slices from the same fragment). Also determines the
        number of reporter slices for each fragment.

        Returns:
         pd.DataFrame: Slices aggregated by parental read name.

        """

        df = (
            self.slices.sort_values(["parent_read", "chrom", "start"])
            .groupby("parent_read", as_index=False, sort=False)
            .agg(
                id=("parent_id", lambda df: df.head(n=1)),
                unique_slices=("slice", "nunique"),
                pe=("pe", "first"),
                mapped=("mapped", "sum"),
                multimapped=("multimapped", "sum"),
                unique_capture_sites=("capture", "nunique"),
                capture_count=("capture_count", "sum"),
                unique_exclusions=("exclusion", "nunique"),
                exclusion_count=("exclusion_count", "sum"),
                unique_restriction_fragments=("restriction_fragment", "nunique"),
                blacklist=("blacklist", "sum"),
                coordinates=("coordinates", "|".join),
            )
        )

        # df["unique_capture_sites"] = (
        #     df["unique_capture_sites"] - 1
        # )  # nunique identifies '.' as a capture site
        # df["unique_exclusions"] = df["unique_exclusions"] - 1  # as above

        # Add the number of reporters to the dataframe.
        # Only consider a reporter if at least one capture slice is present
        # in the fragment.
        df["reporter_count"] = np.where(
            df["capture_count"] > 0,
            df["mapped"]
            - (df["exclusion_count"] + df["capture_count"] + df["blacklist"]),
            0,
        )

        return df

    @property
    def slice_stats(self):
        slices = self.slices.copy()
        if slices.empty:  # Deal with empty dataframe i.e. no valid slices
            for col in slices:
                slices[col] = np.zeros((10,))

        stats_df = slices.agg(
            {
                "slice_name": "nunique",
                "parent_read": "nunique",
                "mapped": "sum",
                "multimapped": "sum",
                "capture": "nunique",
                "capture_count": lambda col: (col > 0).sum(),
                "exclusion_count": lambda col: (col > 0).sum(),
                "blacklist": "sum",
            }
        )

        stats_df = stats_df.rename(
            {
                "slice_name": "unique_slices",
                "parent_read": "unique_fragments",
                "multimapped": "multimapping_slices",
                "capture": "unique_capture_sites",
                "capture_count": "number_of_capture_slices",
                "exclusion_count": "number_of_slices_in_exclusion_region",
                "blacklist": "number_of_slices_in_blacklisted_region",
            }
        )

        return stats_df

    @property
    def frag_stats(self) -> pd.DataFrame:
        """
        Statistics aggregated at the fragment level.

        As this involves slice aggregation it can be rather slow
        for large datasets. It is recomended to only use this
        property if it is required.


        Returns:
         pd.DataFrame: Fragment level statistics
        """

        return self.fragments.agg(
            {
                "parent_read": "nunique",
                "mapped": lambda col: (col > 1).sum(),
                "multimapped": lambda col: (col > 0).sum(),
                "capture_count": lambda col: (col > 0).sum(),
                "exclusion_count": lambda col: (col > 0).sum(),
                "blacklisted_slices": lambda col: (col > 0).sum(),
                "reporter_count": lambda col: (col > 0).sum(),
            }
        ).rename(
            {
                "parent_read": "unique_fragments",
                "multimapped": "fragments_with_multimapping_slices",
                "capture_count": "fragments_with_capture_sites",
                "exclusion_count": "fragments_with_excluded_regions",
                "blacklisted_slices": "fragments_with_blacklisted_regions",
                "reporter_count": "fragments_with_reporter_slices",
            }
        )

    @property
    def reporters(self) -> pd.DataFrame:
        # Return any slice with a  N/A value
        return self.slices.query("capture_count < 1")

    @property
    def captures(self) -> pd.DataFrame:
        """
        Extracts capture slices from slices dataframe

        i.e. slices that do not have a null capture name

        Returns:
         pd.DataFrame: Capture slices

        """
        # Return any slice with a non N/A capture value
        return self.slices.query("capture_count == 1")

    @property
    def capture_site_stats(self) -> pd.Series:
        """Extracts the number of unique capture sites."""
        return self.captures["capture"].value_counts()

    @property
    def merged_captures_and_reporters(self) -> pd.DataFrame:
        """
        Merges captures and reporters sharing the same parental id.

        Capture slices and reporter slices with the same parental read id are
        merged together. The prefixes 'capture' and 'reporter' are used to
        identify slices marked as either captures or reporters.

        Returns:
         pd.DataFrame: Merged capture and reporter slices
        """

        captures = (
            self.captures.set_index("parent_read")
            .add_prefix("capture_")
            .rename(columns={"capture_capture": "capture"})
        )

        reporters = self.reporters.set_index("parent_read").add_prefix("reporter_")

        # Join reporters to captures using the parent read name
        captures_and_reporters = captures.join(reporters).reset_index()

        return captures_and_reporters

    @property
    def cis_or_trans_stats(self) -> pd.DataFrame:
        """
        Extracts reporter cis/trans statistics from slices.

        Returns:
         pd.DataFrame: Reporter cis/trans statistics
        """
        cap_and_rep = self.merged_captures_and_reporters.copy()

        cap_and_rep["cis/trans"] = np.where(
            cap_and_rep["capture_chrom"] == cap_and_rep["reporter_chrom"],
            "cis",
            "trans",
        )

        # Aggregate by capture site for reporting

        return (
            cap_and_rep.groupby(["capture", "cis/trans"])
            .size()
            .reset_index()
            .rename(columns={"capture": "viewpoint", 0: "count"})
            .assign(sample=self.sample_name, read_type=self.read_type)
        )

    def remove_non_reporter_fragments(self):
        """
        Removes the fragment if it has no reporter slices present (Common)

        """

        fragments_partial = self.slices.groupby("parent_id").agg(
            n_capture=("capture_count", "sum"),
            n_mapped=("mapped", "sum"),
            n_blacklist=("blacklist", "sum"),
            n_exclusions=("exclusion_count", "sum"),
        )

        fragments_with_reporters = fragments_partial.query(
            "(n_mapped - n_capture - n_blacklist - n_exclusions) > 0"
        )

        self.slices = (
            self.slices.set_index("parent_id")
            .loc[fragments_with_reporters.index]
            .reset_index()
        )

    def remove_multi_capture_fragments(self):
        """
        Removes double capture fragments.

        All slices (i.e. the entire fragment) are removed if more than
        one capture probe is present i.e. a double capture (V. Common)

        """
        fragments_n_captures = self.slices.groupby("parent_id")["capture"].nunique()
        single_capture_fragments = fragments_n_captures[fragments_n_captures == 1]

        self.slices = (
            self.slices.set_index("parent_id")
            .loc[single_capture_fragments.index]
            .reset_index()
        )

    def remove_viewpoint_adjacent_restriction_fragments(self, n_adjacent: int = 1):
        """
        Deals with an odd situation in which a reporter spanning two adjacent capture sites is not removed.

        Example:
         ------Capture 1----/------Capture 2------\
                  -----REP--------

        In this case the "reporter" slice is not considered either a capture or exclusion.

        These cases are dealt with by explicitly removing reporters on restriction fragments
        adjacent to capture sites.

        Args:
         n_adjacent: Number of adjacent restriction fragments to remove

        """

        slices_with_viewpoint = self.slices_with_viewpoint

        # Create a per viewpoint dataframe of adjacent fragment ranges
        restriction_fragments_viewpoint = (
            self.captures.set_index("capture")["restriction_fragment"]
            .drop_duplicates()
            .reset_index()
            .assign(
                exclusion_start=lambda df: df["restriction_fragment"] - n_adjacent,
                exclusion_end=lambda df: df["restriction_fragment"] + n_adjacent,
            )
        )

        slices_with_viewpoint = slices_with_viewpoint.merge(
            restriction_fragments_viewpoint[
                ["capture", "exclusion_start", "exclusion_end"]
            ],
            left_on="viewpoint",
            right_on="capture",
        )

        # Mark slices between the exclusion zones but ignore capture slices
        excluded_slices = slices_with_viewpoint.query(
            "(exclusion_start <= restriction_fragment <= exclusion_end) and (capture_count == 0)"
        )

        self.slices = (
            self.slices.loc[lambda df: ~df["parent_id"].isin(excluded_slices["parent_id"])]
        )


class TriCSliceFilter(CCSliceFilter):
    """
    Perform Tri-C slice filtering (inplace) and reporter identification.

    SliceFilter tuned specifically for Tri-C data. Whilst the vast majority of filters
    are inherited from CCSliceFilter, this class has addtional methods for Tri-C analysis
    i.e. remove_slices_with_one_reporter. The default filtering order is:

     - remove_unmapped_slices
     - remove_slices_without_re_frag_assigned
     - remove_orphan_slices
     - remove_multi_capture_fragments
     - remove_blacklisted_slices
     - remove_non_reporter_fragments
     - remove_viewpoint_adjacent_restriction_fragments
     - remove_duplicate_re_frags
     - remove_duplicate_slices
     - remove_duplicate_slices_pe
     - remove_non_reporter_fragments
     - remove_slices_with_one_reporter

    See the individual methods for further details.

    Attributes:
     slices (pd.DataFrame): Annotated slices dataframe.
     fragments (pd.DataFrame): Slices dataframe aggregated by parental read.
     reporters (pd.DataFrame): Slices identified as reporters.
     filter_stages (dict): Dictionary containg stages and a list of class methods (str) required to get to this stage.
     slice_stats (pd.DataFrame): Provides slice level statistics.
     read_stats (pd.DataFrame): Provides statistics of slice filtering at the parental read level.
     filter_stats (pd.DataFrame): Provides statistics of read filtering."""

    def __init__(self, slices, filter_stages=None, **sample_kwargs):

        if filter_stages:
            self.filter_stages = filter_stages
        else:
            filter_stages = {
                "pre-filtering": [
                    "get_unfiltered_slices",
                ],
                "mapped": [
                    "remove_unmapped_slices",
                    "remove_slices_without_re_frag_assigned",
                ],
                "contains_single_capture": [
                    "remove_orphan_slices",
                    "remove_multi_capture_fragments",
                ],
                "contains_capture_and_reporter": [
                    "remove_blacklisted_slices",
                    "remove_non_reporter_fragments",
                ],
                "duplicate_filtered": [
                    "remove_duplicate_re_frags",
                    "remove_duplicate_slices",
                    "remove_duplicate_slices_pe",
                    "remove_non_reporter_fragments",
                ],
                "tric_reporter": ["remove_slices_with_one_reporter"],
            }

        super(TriCSliceFilter, self).__init__(slices, filter_stages, **sample_kwargs)

    def remove_slices_with_one_reporter(self):
        """Removes fragments if they do not contain at least two reporters."""
        fragments_triplets = self.fragments.query("reporter_count > 1")
        self.slices = self.slices.loc[
            lambda df: df["parent_read"].isin(fragments_triplets["parent_read"])
        ]


class TiledCSliceFilter(SliceFilter):
    """
    Perform Tiled-C slice filtering (inplace) and reporter identification.

    SliceFilter tuned specifically for Tiled-C data. This class has addtional methods
    to remove common artifacts in Tiled-C data i.e. non-capture fragments,
    multi-capture (with different tiled regions) fragments.
    A reporter is defined differently in a Tiled-C analysis as a reporter slice can also
    be a capture slice.

    The default filter order is as follows:

     - remove_unmapped_slices
     - remove_orphan_slices
     - remove_blacklisted_slices
     - remove_non_capture_fragments
     - remove_dual_capture_fragments
     - remove_slices_without_re_frag_assigned
     - remove_duplicate_re_frags
     - remove_duplicate_slices
     - remove_duplicate_slices_pe
     - remove_orphan_slices

    See the individual methods for further details.

    Attributes:
     slices (pd.DataFrame): Annotated slices dataframe.
     fragments (pd.DataFrame): Slices dataframe aggregated by parental read.
     reporters (pd.DataFrame): Slices identified as reporters.
     filter_stages (dict): Dictionary containg stages and a list of class methods (str) required to get to this stage.
     slice_stats (pd.DataFrame): Provides slice level statistics.
     read_stats (pd.DataFrame): Provides statistics of slice filtering at the parental read level.
     filter_stats (pd.DataFrame): Provides statistics of read filtering.

    """

    def __init__(self, slices, filter_stages=None, **sample_kwargs):

        if not filter_stages:
            filter_stages = {
                "pre-filtering": [
                    "get_unfiltered_slices",
                ],
                "mapped": ["remove_unmapped_slices", "remove_orphan_slices"],
                "not_blacklisted": ["remove_blacklisted_slices"],
                "contains_capture": [
                    "remove_non_capture_fragments",
                    "remove_dual_capture_fragments",
                ],
                "duplicate_filtered": [
                    "remove_slices_without_re_frag_assigned",
                    "remove_duplicate_re_frags",
                    "remove_duplicate_slices",
                    "remove_duplicate_slices_pe",
                ],
                "has_reporter": ["remove_orphan_slices"],
            }

        super(TiledCSliceFilter, self).__init__(slices, filter_stages, **sample_kwargs)

    @property
    def fragments(self) -> pd.DataFrame:
        df = (
            self.slices.sort_values(["parent_read", "chrom", "start"])
            .groupby("parent_read", as_index=False, sort=False)
            .agg(
                id=("parent_id", "first"),
                unique_slices=("slice", "nunique"),
                pe=("pe", "first"),
                mapped=("mapped", "sum"),
                multimapped=("multimapped", "sum"),
                capture_count=("capture_count", "sum"),
                unique_restriction_fragments=("restriction_fragment", "nunique"),
                blacklisted_slices=("blacklist", "sum"),
                coordinates=("coordinates", "|".join),
            )
        )

        return df

    @property
    def slice_stats(self):
        stats_df = self.slices.agg(
            {
                "slice_name": "nunique",
                "parent_read": "nunique",
                "mapped": "sum",
                "multimapped": "sum",
                "capture_count": lambda col: (col > 0).sum(),
                "blacklist": "sum",
            }
        )

        stats_df = stats_df.rename(
            {
                "slice_name": "unique_slices",
                "parent_read": "unique_fragments",
                "multimapped": "multimapping_slices",
                "capture_count": "number_of_capture_slices",
                "blacklist": "number_of_slices_in_blacklisted_region",
            }
        )

        return stats_df

    @property
    def cis_or_trans_stats(self) -> pd.DataFrame:
        """
        Extracts reporter cis/trans statistics from slices.

        Unlike Capture-C/Tri-C reporter slice can also be capture slices as
        all slices within the capture region are considered as reporters. To extract
        cis/trans statistics, one capture slice in each fragment is considered to be
        the "primary capture" this then enables merging of this "primary capture" with
        the other reporters both inside and outside of the tiled region.

        Returns:
         pd.DataFrame: Reporter cis/trans statistics
        """

        interactions_by_capture = dict()

        for capture_site, df_cap in self.slices.query("capture_count == 1").groupby(
            "capture"
        ):

            capture_chrom = df_cap.iloc[0]["chrom"]
            df_primary_capture = df_cap.groupby(
                "parent_read"
            ).first()  # Artifact required as need to call one slice the "capture"
            df_not_primary_capture = df_cap.loc[
                ~(df_cap["slice_name"].isin(df_primary_capture["slice_name"]))
            ]
            df_outside_capture = self.slices.query("capture != capture").loc[
                lambda df_rep: df_rep["parent_read"].isin(df_cap["parent_read"])
            ]

            df_pseudo_reporters = pd.concat(
                [df_not_primary_capture, df_outside_capture]
            )
            n_cis_interactions = df_pseudo_reporters.query(
                f'chrom == "{capture_chrom}"'
            ).shape[0]
            n_trans_interactions = df_pseudo_reporters.shape[0] - n_cis_interactions

            interactions_by_capture[capture_site] = {
                "cis": n_cis_interactions,
                "trans": n_trans_interactions,
            }

        return (
            pd.DataFrame(interactions_by_capture)
            .transpose()
            .reset_index()
            .rename(columns={"index": "capture"})
            .melt(id_vars="capture", var_name="cis/trans", value_name="count")
            .sort_values("capture")
            .assign(sample=self.sample_name, read_type=self.read_type)
            .rename(columns={"capture": "viewpoint"})
        )

    def remove_slices_outside_capture(self):
        """Removes slices outside of capture region(s)"""
        self.slices = self.slices.query("capture_count != 0")

    def remove_non_capture_fragments(self):
        """Removes fragments without a capture assigned"""
        fragments_with_capture = (
            self.slices.groupby("parent_id")["capture_count"]
            .sum()
            .reset_index()
            .query("capture_count > 0")
        )
        self.slices = self.slices[
            self.slices["parent_read"].isin(fragments_with_capture["parent_id"])
        ]

    def remove_dual_capture_fragments(self):
        """
        Removes a fragment with multiple different capture sites.

        Modified for TiledC filtering as the fragment dataframe is generated
        slightly differently.
        """
        multicapture_fragments = (
            self.slices.query("capture_count == 1")
            .groupby("parent_read")["capture"]
            .nunique()
            > 1
        )
        self.slices = (
            self.slices.set_index("parent_read")
            .loc[~multicapture_fragments]
            .reset_index()
        )
