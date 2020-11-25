import os
import time
from datetime import timedelta
from functools import wraps

import numpy as np
import pandas as pd
import pysam


def get_timing(task_name=None):
    """Decorator:
    Gets the time taken by the wrapped function
    """

    def wrapper(f):
        @wraps(f)
        def wrapped(*args, **kwargs):
            time_start = time.perf_counter()
            result = f(*args, **kwargs)
            time_end = time.perf_counter()

            time_taken = timedelta(seconds=(time_end - time_start))
            print(f"Completed {task_name} in {time_taken} (hh:mm:ss.ms)")
            return result

        return wrapped

    return wrapper


def parse_alignment(aln):
    """Parses reads from a bam file into a list.

    Args:
      aln: pysam.AlignmentFile
    Returns:
      List containing:
      -read name
      -parent reads
      -flashed status
      -slice number
      -mapped status
      -multimapping status
      -chromosome number (e.g. chr10)
      -start (e.g. 1000)
      -end (e.g. 2000)
      -coords e.g. (chr10:1000-2000)
    """

    read_name = aln.query_name
    parent_read, pe, slice_number, uid = read_name.split("|")
    ref_name = aln.reference_name
    ref_start = aln.reference_start
    ref_end = aln.reference_end
    # Check if read mapped
    if aln.is_unmapped:
        mapped = 0
        multimapped = 0
        ref_name = "unmapped"
        ref_start = ""
        ref_end = ""
        coords = ""
    else:
        mapped = 1
        coords = f"{ref_name}:{ref_start}-{ref_end}"
        # Check if multimapped
        if aln.is_secondary:
            multimapped = 1
        else:
            multimapped = 0
    return [
        read_name,
        parent_read,
        pe,
        slice_number,
        uid,
        mapped,
        multimapped,
        ref_name,
        ref_start,
        ref_end,
        coords,
    ]


@get_timing(task_name="processing BAM file")
def parse_bam(bam):
    """Uses parse_alignment function convert bam file to a dataframe.

    Args:
     bam: File name of bam file to process.

    Returns:
     Dataframe with columns:
     -'read_name'
     -'parent_read'
     -'pe'
     -'slice'
     -'mapped'
     -'multimapped'
     -'chrom'
     -'start'
     -'end'
     -'coordinates'
    """

    df_bam = pd.DataFrame(
        [
            parse_alignment(aln)
            for aln in pysam.AlignmentFile(bam, "rb").fetch(until_eof=True)
        ],
        columns=[
            "read_name",
            "parent_read",
            "pe",
            "slice",
            "uid",
            "mapped",
            "multimapped",
            "chrom",
            "start",
            "end",
            "coordinates",
        ],
    )
    df_bam.set_index(["read_name", "chrom", "start"], inplace=True)
    return df_bam


@get_timing(task_name="merging annotations with BAM input")
def merge_annotations(df, annotations):
    """Combines annotations with the parsed bam file output.

    Uses pandas outer join on the indexes to merge annotations
    e.g. number of capture probe overlaps.

    Annotation tsv must have the index as the first column and this index
    must have intersecting keys with the first dataframe's index.

    Args:
     df: pd.Dataframe to merge with annotations
     annotations: Filename of .tsv to read and merge with df

    Returns:
     Merged dataframe

    """

    # Update, now using chrom and start to stop issues with multimapping
    df_ann = pd.read_csv(
        annotations, sep="\t", header=0, index_col=["read_name", "chrom", "start"]
    )
    df_ann = df_ann.drop(columns="end", errors="ignore")

    return (
        df.join(df_ann, how="inner")
        .drop(columns=["read_name.1"], errors="ignore")
        .reset_index()
    )


class SliceFilter:
    def __init__(self, slices, filter_stages=None):

        self.slices = slices.copy()

        if filter_stages:
            self.filter_stages = filter_stages
        else:
            raise ValueError("Filter stages not provided")

        self.filtered = False
        self.filter_stats = pd.DataFrame()

    @property
    def slice_stats(self):
        raise NotImplementedError("Override this method")

    @property
    def fragments(self):
        df = (
            self.slices.sort_values(["parent_read", "chrom", "start"])
            .groupby("parent_read", as_index=False)
            .agg(
                {
                    "slice": "nunique",
                    "pe": "first",
                    "mapped": "sum",
                    "multimapped": "sum",
                    "exclusion": "nunique",
                    "exclusion_count": "sum",
                    "restriction_fragment": "nunique",
                    "blacklist": "sum",
                    "coordinates": "|".join,
                }
            )
        )

    @property
    def reporters(self):
        raise NotImplementedError("Override this property")

    def filter_slices(self, output_slices=False, output_location="."):

        for stage, filters in self.filter_stages.items():
            for filt in filters:
                # Call all of the filters in the filter_stages dict in order
                print(f"Filtering: {filt}")
                getattr(self, filt)()  # Gets and calls the selected method
                print(f"Number of slices: {self.slices.shape[0]}")

                if output_slices == "filter":
                    self.slices.to_csv(os.path.join(output_location, f"{filt}.tsv.gz"))

            if output_slices == "stage":
                self.slices.to_csv(os.path.join(output_location, f"{stage}.tsv.gz"))

            self.filter_stats[
                stage
            ] = self.slice_stats  # Store the stats for each stage

    def modify_re_frag(self, frag: str, adjust=1):
        """Increases/Decreases the RE frag number.

        e.g. modify_re_frag(DpnII_chr10_5, adjust=1) -> DpnII_chr10_6

        Args:
         frag: Name of restriction fragment (str)
         adjust: Adjust fragment identifier number by value

        Returns:
         Modified fragment name (str)
        """
        if frag != ".":
            enzyme, chrom, index = frag.split("_")
            return "_".join([enzyme, chrom, str(int(index) + adjust)])

    def remove_unmapped_slices(self):
        """Removes slices marked as unmapped (Uncommon)

        Returns:
         CCSliceFilter
        """
        self.slices = self.slices.query("mapped == 1")

    def remove_orphan_slices(self):
        """Remove fragments with only one aligned slice (Common)

        Returns:
         CCSliceFilter
        """
        # TODO: Modify this function to speed it up (remove groupby and replace with ['parent_read].duplicated())

        fragments = self.fragments
        fragments_multislice = fragments.query("unique_slices > 1")
        self.slices = self.slices[
            self.slices["parent_read"].isin(fragments_multislice["parent_read"])
        ]

        # TODO: Fix this
        #self.slices = self.slices.loc[self.slices['parent_read'].duplicated()]


    def remove_duplicate_re_frags(self):
        """Prevent the same restriction fragment being counted more than once (Uncommon).
        i.e. --RE_FRAG1--\----Capture----\---RE_FRAG1----

        Returns:
          CCSliceFilter

        """
        self.slices = self.slices.sample(frac=1).drop_duplicates(
            subset=["parent_read", "restriction_fragment"], keep="first"
        )

    def remove_slices_without_re_frag_assigned(self):
        self.slices = self.slices.query('restriction_fragment != "."')

    def remove_duplicate_slices(self):
        """Remove all slices if the slice coordinates and slice order are shared
        with another fragment i.e. are PCR duplicates (Common).

        e.g
                          coordinates
        | Frag 1:  chr1:1000-1250 chr1:1500-1750
        | Frag 2:  chr1:1000-1250 chr1:1500-1750
        | Frag 3:  chr1:1050-1275 chr1:1600-1755
        | Frag 4:  chr1:1500-1750 chr1:1000-1250

        Frag 2 removed. Frag 1,3,4 retained

        Returns:
         CCSliceFilter
        """

        frags_deduplicated = self.fragments.sample(frac=1).drop_duplicates(
            subset="coordinates", keep="first"
        )

        self.slices = self.slices[
            self.slices["parent_read"].isin(frags_deduplicated["parent_read"])
        ]

    def remove_duplicate_slices_pe(self):
        """Removes PCR duplicates from non-flashed (PE) fragments (Common).
        Sequence quality is often lower at the 3' end of reads leading to variance in mapping coordinates.
        PCR duplicates are removed by checking that the fragment start and end are not duplicated in the dataframe.

        Returns:
         CCSliceFilter
        """
        if (
            self.slices["pe"].str.contains("unflashed").sum() > 1
        ):  # at least one un-flashed
            fragments = self.fragments.assign(
                read_start=lambda df: df["coordinates"]
                .str.split("|")
                .str[0]
                .str.split(r":|-")
                .str[1],
                read_end=lambda df: df["coordinates"]
                .str.split("|")
                .str[-1]
                .str.split(r":|-")
                .str[-1],
            )

            fragments_pe = fragments.query('pe == "unflashed"')
            fragments_pe_duplicated = fragments_pe[
                fragments_pe.duplicated(subset=["read_start", "read_end"])
            ]

            self.slices = self.slices[
                ~(
                    self.slices["parent_read"].isin(
                        fragments_pe_duplicated["parent_read"]
                    )
                )
            ]  # Slices not in duplicated

    def remove_excluded_slices(self):
        """Removes any slices in the exclusion region (default 1kb) and a blacklist (if supplied) (V. Common)

        Returns:
         CCSliceFilter"""
        self.slices = self.slices.query("exclusion_count < 1")

    def remove_blacklisted_slices(self):
        self.slices = self.slices.query("blacklist < 1")


class CCSliceFilter(SliceFilter):
    """Class containing methods for filtering slices and reporting
    slice/fragment statistics.

    Attributes:
         slices: DataFrame containing aligned reads and annotations.
         filter_stages: Dictionary containing name of stage (for the purpose of storing statistics e.g. duplicate filtered)
                        together with the additional filter names (present in the class) to reach this stage.

    Slices DataFrame must have the following columns:

    - read_name: Unique aligned read identifier (e.g. XZKG:889:11|flashed|1)
    - parent_read: Identifier shared by slices from same fragment (e.g.XZKG:889:11)
    - pe: Read combined by FLASh or not (i.e. "flashed" or "pe")
    - mapped: Alignment is mapped (e.g. 0/1)
    - slice: Slice number (e.g. 0)
    - capture: Capture site intersecting slice (e.g. Slc25A37)
    - capture_count: Number of capture probes overlapping slice (e.g. 1)
    - exclusion: Read present in excluded region (e.g. Slc25A37)
    - exclusion_count: Number of excluded regions overlapping slice (e.g. 1)
    - blacklist: Read present in excluded region (e.g. 0)
    - coordinates: Genome coordinates (e.g. chr1|1000|2000)

    """

    def __init__(self, slices, filter_stages=None):

        if not filter_stages:
            filter_stages = {
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
                    "remove_multicapture_reporters",
                ],
                "duplicate_filtered": [
                    "remove_slices_without_re_frag_assigned",
                    "remove_duplicate_re_frags",
                    "remove_duplicate_slices",
                    "remove_duplicate_slices_pe",
                    "remove_non_reporter_fragments",
                ],
            }

        super(CCSliceFilter, self).__init__(slices, filter_stages)

    @property
    def fragments(self):
        """Summarises slices at the fragment level.

         Uses pandas groupby to aggregate slices by their parental read name
         (shared by all slices from the same fragment). Also determines the
         number of reporter slices for each fragment.

        Returns:
          Dataframe of slices aggregated by fragment

        """
        df = (
            self.slices.sort_values(["parent_read", "chrom", "start"])
            .groupby("parent_read", as_index=False)
            .agg(
                {
                    "slice": "nunique",
                    "pe": "first",
                    "mapped": "sum",
                    "multimapped": "sum",
                    "capture": "nunique",
                    "capture_count": "sum",
                    "exclusion": "nunique",
                    "exclusion_count": "sum",
                    "restriction_fragment": "nunique",
                    "blacklist": "sum",
                    "coordinates": "|".join,
                }
            )
        )
        df["capture"] = df["capture"] - 1  # nunique identifies '.' as a capture site
        df["exclusion"] = df["exclusion"] - 1  # as above

        # Add the number of reporters to the dataframe.
        # Only consider a reporter if at least one capture slice is present
        # in the fragment.
        df["reporter_count"] = np.where(
            df["capture_count"] > 0,
            df["mapped"]
            - (df["exclusion_count"] + df["capture_count"] + df["blacklist"]),
            0,
        )

        # Rename for clarity
        df = df.rename(
            columns={
                "capture": "unique_capture_sites",
                "exclusion": "unique_exclusion_sites",
                "restriction_fragment": "unique_restriction_fragments",
                "slice": "unique_slices",
                "blacklist": "blacklisted_slices",
            }
        )
        return df

    @property
    def slice_stats(self):
        """Gets statisics at a slice level.

        Aggregates slices to determine the number of:
        -unique slices
        -unique fragments
        -unique capture sites
        -capture slices
        -excluded slices
        -blacklisted slices

        Returns:
         Dataframe containing slice statistics
        """
        stats_df = self.slices.agg(
            {
                "read_name": "nunique",
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
                "read_name": "unique_slices",
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
    def frag_stats(self):
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
    def reporters(self):
        """Extracts reporter slices from slices dataframe
        i.e. non-capture slices

        Returns:
         Dataframe containg all non-capture slices"""
        return self.slices.query('capture == "."')

    @property
    def captures(self):
        """Extracts capture slices from slices dataframe
        i.e. slices that do not have a null capture name

        Returns:
         Dataframe containg all capture slices"""
        return self.slices.query('~(capture == ".")')

    @property
    def capture_site_stats(self):
        return self.captures["capture"].value_counts()

    @property
    def merged_captures_and_reporters(self):

        captures = (
            self.captures.set_index("parent_read")
            .add_prefix("capture_")
            .rename(columns={"capture_capture": "capture"})
        )

        reporters = self.reporters.set_index("parent_read").add_prefix("reporter_")

        # Join reporters to captures using the parent read name
        captures_and_reporters = (
            captures.join(reporters).dropna(axis=0, how="any").reset_index()
        )

        return captures_and_reporters

    @property
    def cis_or_trans_stats(self):
        cap_and_rep = self.merged_captures_and_reporters.copy()

        cap_and_rep["cis/trans"] = np.where(
            cap_and_rep["capture_chrom"] == cap_and_rep["reporter_chrom"],
            "cis",
            "trans",
        )

        # Aggregate by capture site for reporting
        interactions_by_capture = pd.DataFrame(
            cap_and_rep.groupby("capture")["cis/trans"].value_counts()
        )

        return interactions_by_capture

    def remove_non_reporter_fragments(self):
        """Removes all slices (i.e. the entire fragment) if it has no reporter slices present (Common)

        Returns:
         CCSliceFilter
        """
        frags_reporter = self.fragments.query("reporter_count > 0")
        self.slices = self.slices[
            self.slices["parent_read"].isin(frags_reporter["parent_read"])
        ]

    def remove_multi_capture_fragments(self):
        """Removes all slices (i.e. the entire fragment) if more than
        one capture probe is present i.e. double captures (V. Common)

        Returns:
         CCSliceFilter
        """
        frags_capture = self.fragments.query("0 < unique_capture_sites < 2")
        self.slices = self.slices[
            self.slices["parent_read"].isin(frags_capture["parent_read"])
        ]

    def remove_multicapture_reporters(self, n_adjacent=1):
        """| Deals with an odd situation in which a reporter spanning two adjacent capture sites is not removed.
        | e.g.
        | ------Capture 1----/------Capture 2------
        |                      -----REP--------
        |
        | In this case the "reporter" slice is not considered either a capture or exclusion.

        | These cases are dealt with by explicitly removing reporters on restriction fragments
        | adjacent to capture sites.

        | The number of adjacent RE fragments can be adjusted with n_adjacent.

        | Returns:
        |  CCSliceFilter
        """

        captures = self.captures
        re_frags = captures["restriction_fragment"].unique()

        # Generates a list of restriction fragments to be excluded from further analysis
        excluded_fragments = [
            self.modify_re_frag(frag, modifier)
            for frag in re_frags
            for modifier in range(-n_adjacent, n_adjacent + 1)
        ]

        # Remove non-capture slices (reporters) in excluded regions
        self.slices = self.slices[
            (self.slices["capture_count"] > 0)
            | (~self.slices["restriction_fragment"].isin(excluded_fragments))
        ]


class TriCSliceFilter(CCSliceFilter):
    def __init__(self, slices, filter_stages=None):

        if filter_stages:
            self.filter_stages = filter_stages
        else:
            filter_stages = {
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

        super(TriCSliceFilter, self).__init__(slices, filter_stages)

    def remove_slices_with_one_reporter(self):
        fragments_triplets = self.fragments.query("reporter_count > 1")
        self.slices = self.slices.loc[
            lambda df: df["parent_read"].isin(fragments_triplets["parent_read"])
        ]


class TiledCSliceFilter(SliceFilter):
    def __init__(self, slices, filter_stages=None):

        if not filter_stages:
            filter_stages = {
                "mapped": ["remove_unmapped_slices", "remove_orphan_slices"],
                "not_blacklisted": ["remove_blacklisted_slices"],
                "contains_capture": ["remove_non_capture_fragments"],
                "duplicate_filtered": [
                    "remove_slices_without_re_frag_assigned",
                    "remove_duplicate_re_frags",
                    "remove_duplicate_slices",
                    "remove_duplicate_slices_pe",
                ],
                "has_reporter": ["remove_orphan_slices"],
            }

        super(TiledCSliceFilter, self).__init__(slices, filter_stages)

    @property
    def fragments(self):
        """Summarises slices at the fragment level.

         Uses pandas groupby to aggregate slices by their parental read name
         (shared by all slices from the same fragment). Also determines the
         number of reporter slices for each fragment.

        Returns:
          Dataframe of slices aggregated by fragment

        """
        df = (
            self.slices.sort_values(["parent_read", "chrom", "start"])
            .groupby("parent_read", as_index=False)
            .agg(
                {
                    "slice": "nunique",
                    "pe": "first",
                    "mapped": "sum",
                    "multimapped": "sum",
                    "capture_count": "sum",
                    "restriction_fragment": "nunique",
                    "blacklist": "sum",
                    "coordinates": "|".join,
                }
            )
        )

        # Rename for clarity
        df = df.rename(
            columns={
                "restriction_fragment": "unique_restriction_fragments",
                "slice": "unique_slices",
                "blacklist": "blacklisted_slices",
            }
        )
        return df

    @property
    def slice_stats(self):
        """Gets statisics at a slice level.

        Aggregates slices to determine the number of:
        -unique slices
        -unique fragments
        -unique capture sites
        -capture slices
        -excluded slices
        -blacklisted slices

        Returns:
         Dataframe containing slice statistics
        """
        stats_df = self.slices.agg(
            {
                "read_name": "nunique",
                "parent_read": "nunique",
                "mapped": "sum",
                "multimapped": "sum",
                "capture_count": lambda col: (col > 0).sum(),
                "blacklist": "sum",
            }
        )

        stats_df = stats_df.rename(
            {
                "read_name": "unique_slices",
                "parent_read": "unique_fragments",
                "multimapped": "multimapping_slices",
                "capture_count": "number_of_capture_slices",
                "blacklist": "number_of_slices_in_blacklisted_region",
            }
        )

        return stats_df

    def remove_slices_outside_capture(self):
        self.slices = self.slices.query('capture != "."')

    def remove_non_capture_fragments(self):
        fragments_with_capture = self.fragments.query("capture_count > 0")
        self.slices = self.slices[
            self.slices["parent_read"].isin(fragments_with_capture["parent_read"])
        ]


@get_timing(task_name="analysis of bam file")
def main(input_bam, annotations, output_prefix, stats_output, method="capture"):

    # TODO: Make sure that the dataframe is correct before processing check all the required columns are present. Check that '.' is the default na value for strings

    # Read bam file and merege annotations
    df_alignment = parse_bam(input_bam)
    df_alignment = merge_annotations(df_alignment, annotations)

    slice_filters_dict = {"capture": CCSliceFilter, "tri": TriCSliceFilter}

    if method == "capture" or method == "tri":
        # Initialise SliceFilter with default args
        slice_filter = slice_filters_dict[method](slices=df_alignment)

        # Filter slices using the slice_filter
        slice_filter.filter_slices()

        # Save filtering statisics
        slice_filter.filter_stats.to_csv(f"{stats_output}.slice.stats", sep="\t")

        # Save reporter stats
        slice_filter.cis_or_trans_stats.to_csv(
            f"{stats_output}.reporter.stats", sep="\t"
        )

        # Output fragments
        slice_filter.fragments.to_csv(
            f"{output_prefix}.fragments.tsv.gz", sep="\t", index=False
        )

        # Output the reporter DataFrame for each capture site
        for capture_site, df_rep in slice_filter.merged_captures_and_reporters.groupby(
            "capture"
        ):
            df_rep.to_csv(
                f"{ output_prefix}.{capture_site}.tsv.gz", sep="\t", index=False
            )

    elif method == "tiled":

        # Initialise TiledCSliceFilter class
        slice_filter = TiledCSliceFilter(df_alignment)

        # Filter slices TiledC-style
        slice_filter.filter_slices()

        # Output fragments
        slice_filter.fragments.to_csv(
            f"{output_prefix}.fragments.tsv.gz", index=None, sep="\t"
        )

        # Save filtering statisics
        slice_filter.filter_stats.to_csv(f"{stats_output}.slice.stats", sep="\t")

        # Output filtered slices
        reporters = slice_filter.slices.add_prefix("reporter_").rename(
            columns={
                "reporter_parent_read": "parent_read",
                "reporter_capture": "capture",
            }
        )

        # Assumes only one region tiled per assay
        tiled_region_name = reporters.query('capture != "."')['capture'].unique()[0]
        reporters.to_csv(
            f"{ output_prefix}.{tiled_region_name}.tsv.gz", sep="\t", index=False
        )

        # Output dummy reporter stats
        df_reporter_stats = pd.DataFrame(
            {
                "capture": [
                    tiled_region_name,
                ],
                "cis/trans": [
                    "cis",
                ],
                "count": [
                    slice_filter.filter_stats.loc[
                        "number_of_capture_slices", "has_reporter"
                    ],
                ],
            }
        )

        df_reporter_stats.to_csv(
            f"{stats_output}.reporter.stats", sep="\t", index=False
        )
