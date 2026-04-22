import warnings
from typing import Union, List, Literal

import pandas as pd
import pyranges1 as pr
import numpy as np
from pandas.api.types import is_categorical_dtype, is_numeric_dtype

from capcruncher.utils import convert_bed_to_pr


def _pyranges_frame(gr: pr.PyRanges) -> pd.DataFrame:
    return pd.DataFrame(gr).copy()

warnings.simplefilter("ignore", category=RuntimeWarning)


def increase_cis_slice_priority(df: pd.DataFrame, score_multiplier: float = 2):
    """
    Prioritizes cis slices by increasing the mapping score.
    """

    df["parent_name"] = df["name"].str.split("|").str[0]

    df_chrom_counts = (
        df[["parent_name", "chrom"]].value_counts().to_frame("slices_per_chrom")
    )
    modal_chrom = (
        df_chrom_counts.groupby("parent_name")["slices_per_chrom"]
        .transform("max")
        .reset_index()
        .set_index("parent_name")["chrom"]
        .to_dict()
    )
    df["fragment_chrom"] = df["parent_name"].map(modal_chrom)
    df["score"] = np.where(
        df["chrom"] == df["fragment_chrom"],
        df["score"] * score_multiplier,
        df["score"] / score_multiplier,
    )

    return df.drop(columns="parent_name")


def remove_duplicates_from_bed(
    bed: pr.PyRanges,
    prioritize_cis_slices: bool = False,
    chroms_to_prioritize: Union[list, np.ndarray] = None,
) -> pr.PyRanges:
    """
    Removes duplicate entries from a PyRanges object.

    Args:
        bed (pr.PyRanges): PyRanges object to be deduplicated.
        prioritize_cis_slices (bool, optional): Prioritize cis slices by increasing the mapping score. Defaults to False.
        chroms_to_prioritize (Union[list, np.ndarray], optional): Chromosomes to prioritize. Defaults to None.

    Returns:
        pr.PyRanges: Deduplicated PyRanges object.
    """

    df = _pyranges_frame(bed).rename(columns=lambda col: col.lower()).rename(
        columns={"chromosome": "chrom"}
    )

    # Shuffle the dataframe to randomize the duplicate removal
    df = df.sample(frac=1)

    if prioritize_cis_slices:
        df = increase_cis_slice_priority(df)

    if "score" in df.columns:
        df = df.sort_values(["score"], ascending=False)

    if chroms_to_prioritize:
        df["is_chrom_priority"] = df["chrom"].isin(chroms_to_prioritize).astype(int)
        df = df.sort_values(["score", "is_chrom_priority"], ascending=False).drop(
            columns="is_chrom_priority"
        )

    return (
        df.drop_duplicates(subset="name", keep="first")
        .sort_values(["chrom", "start"])[["chrom", "start", "end", "name"]]
        .rename(columns=lambda col: col.capitalize())
        .rename(columns={"Chrom": "Chromosome"})
        .pipe(pr.PyRanges)
    )


class Intersection:
    def __init__(
        self,
        bed_a: pr.PyRanges,
        bed_b: pr.PyRanges,
        name: str,
        fraction: float = 0,
        n_cores: int = 1,
    ):
        self.a = bed_a
        self.b = bed_b
        self.name = name
        self.fraction = fraction
        self.n_cores = n_cores

    @property
    def intersection(self) -> pr.PyRanges:
        raise NotImplementedError("Must be implemented in subclass")


class IntersectionGet(Intersection):
    def get_new_dtype(self):
        bed_b = _pyranges_frame(self.b)

        # Determine the dtype of the name column
        if is_numeric_dtype(bed_b["Name"]):
            dtype_new = bed_b["Name"].dtype

        elif is_categorical_dtype(bed_b["Name"]):
            if is_numeric_dtype(bed_b["Name"].cat.categories):
                dtype_new = bed_b["Name"].cat.categories.dtype
                numeric_dtype_mapping = {
                    "int64": "Int64",
                    "float64": "Float64",
                    "float32": "Float32",
                    "int32": "Int32",
                }
                dtype_new = numeric_dtype_mapping.get(str(dtype_new), "Int64")

            else:
                dtype_new = bed_b["Name"].dtype

        else:
            dtype_new = pd.CategoricalDtype([*bed_b["Name"].unique().astype(str)])

        return dtype_new

    @property
    def intersection(self) -> pr.PyRanges:
        dtype_new = self.get_new_dtype()

        df_both = pd.DataFrame(
            self.a.join_overlaps(
                self.b,
                join_type="left",
                report_overlap_column="Overlap",
            )
        ).sort_values("Name")

        # Filter out the non-overlapping regions
        df_both["frac"] = df_both["Overlap"].fillna(0).div(
            df_both["End"] - df_both["Start"]
        )
        df_both[self.name] = np.where(
            df_both["frac"] >= self.fraction, df_both["Name_b"], pd.NA
        )
        df_both[self.name] = df_both[self.name].astype(dtype_new)

        df_both.drop(
            columns=[
                "frac",
                "Overlap",
                "Name_b",
                "Start_b",
                "End_b",
                "Strand_b",
                "Score_b",
            ],
            errors="ignore",
            inplace=True,
        )

        return df_both.pipe(pr.PyRanges)


class IntersectionCount(Intersection):
    @property
    def intersection(self) -> pr.PyRanges:
        joined = pd.DataFrame(
            self.a.join_overlaps(self.b, report_overlap_column="Overlap")
        )

        bed_a = _pyranges_frame(self.a)
        if joined.empty:
            counts = pd.Series(0, index=bed_a["Name"], dtype="int64")
        else:
            joined["fraction"] = joined["Overlap"] / (joined["End"] - joined["Start"])
            counts = (
                joined.loc[joined["fraction"] >= self.fraction]
                .groupby("Name")
                .size()
                .reindex(bed_a["Name"], fill_value=0)
            )

        return pr.PyRanges(
            bed_a.assign(**{self.name: counts.to_numpy(dtype="int64")})
        )


class IntersectionFailed(Intersection):
    @property
    def intersection(self):
        return (
            _pyranges_frame(self.a).assign(**{self.name: pd.NA})
            .assign(**{self.name: lambda df: df[self.name].astype(pd.StringDtype())})
            .pipe(pr.PyRanges)
        )


class BedIntersector:
    def __init__(
        self,
        bed_a: Union[str, pr.PyRanges],
        bed_b: Union[str, pr.PyRanges],
        name: str,
        fraction: float = 0,
        max_cores: int = 1,
    ):
        self.annotation_columns = None

        if isinstance(bed_a, pr.PyRanges):
            self.a = self.process_bed(bed_a)
        elif isinstance(bed_a, (str, pd.DataFrame)):
            self.a = convert_bed_to_pr(bed_a)
            self.a = self.process_bed(self.a)
        else:
            raise ValueError(
                f"bed_a must be of type str, pd.DataFrame, or pr.PyRanges. Got {type(bed_a)}"
            )

        self.b = bed_b if isinstance(bed_b, pr.PyRanges) else convert_bed_to_pr(bed_b)
        self.name = name
        self.fraction = fraction
        self.n_cores = max_cores if len(self.b) > 50_000 else 1

    def get_intersection(self, method: Literal["get", "count"] = "get") -> pr.PyRanges:
        try:
            if self.b.empty:
                _intersection = IntersectionFailed(
                    self.a, self.b, self.name, self.fraction, self.n_cores
                ).intersection
            elif method == "get":
                _intersection = IntersectionGet(
                    self.a, self.b, self.name, self.fraction, self.n_cores
                ).intersection
            elif method == "count":
                _intersection = IntersectionCount(
                    self.a, self.b, self.name, self.fraction, self.n_cores
                ).intersection
            else:
                _intersection = IntersectionFailed(
                    self.a, self.b, self.name, self.fraction, self.n_cores
                ).intersection

        except (
            OSError,
            IndexError,
            FileNotFoundError,
            StopIteration,
            AssertionError,
            ValueError,
        ):
            _intersection = IntersectionFailed(
                self.a, self.b, self.name, self.fraction, self.n_cores
            ).intersection

        # If there are annotation columns, join them to the intersection
        if not self.annotation.empty:
            _intersection = (
                _pyranges_frame(_intersection).set_index("Name")
                .join(self.annotation, how="left")
                .reset_index()
                .pipe(pr.PyRanges)
            )

        # Put the original name back
        _intersection = _pyranges_frame(_intersection).assign(
            Name=lambda df: df.Name.map(self.original_name_mapping)
        )

        return pr.PyRanges(_intersection)

    def process_bed(self, bed: pr.PyRanges):
        # Convert to dataframe
        bed = _pyranges_frame(bed)

        # Create a unique identifier for each slice
        self.uid = pd.util.hash_pandas_object(
            bed.loc[:, ["Chromosome", "Start", "End", "Name"]]
        )
        self.original_name_mapping = dict(zip(self.uid, bed.Name))

        # Add the unique identifier to the bed
        bed = bed.assign(Name=self.uid)

        # Identify colunms that have annotation information
        self.annotation_col_names = [
            col
            for col in bed.columns
            if col not in ["Chromosome", "Start", "End", "Strand", "Score", "Name"]
        ]

        # If there are annotation columns, store them in a separate dataframe
        if self.annotation_col_names:
            self.annotation = bed.set_index("Name").loc[:, self.annotation_col_names]
        else:
            self.annotation = pd.DataFrame()

        # Re-generate the pyranges object with the unique identifier
        bed = bed.loc[:, ["Chromosome", "Start", "End", "Name"]]

        return bed.pipe(pr.PyRanges)


# @ray.remote
# class BedFileIntersection:
#     """
#     Intersect two bed files and return the intersection as a pandas series.

#     Args:
#         bed_a (Union[str, pybedtools.BedTool, pr.PyRanges]): First bed file to intersect.
#         bed_b (Union[str, pybedtools.BedTool, pr.PyRanges]): Second bed file to intersect.
#         name (str, optional): Name of the intersection. Defaults to "b".
#         action (str, optional): Method to use for intersection. Defaults to "get".
#         fraction (float, optional): Minimum fraction of overlap to consider a hit. Defaults to 1e-9.
#     """

#     def __init__(
#         self,
#         bed_a: Union[str, pybedtools.BedTool, pr.PyRanges],
#         bed_b: Union[str, pybedtools.BedTool, pr.PyRanges],
#         name: str = "b",
#         action: str = "get",
#         fraction: float = 1e-9,
#     ):

#         self.a = bed_a
#         self.b = bed_b
#         self.name = name
#         self.action = action
#         self.fraction = fraction

#         self.pr_a = convert_bed_to_pr(self.a, ignore_ray_objrefs=True)

#         import logging

#         logging.basicConfig(level=logging.INFO)
#         self.logger = logging.getLogger(__name__)

#     def _get_intersection(self, pr_b: pr.PyRanges):

#         intersection = (
#             self.pr_a.join(
#                 pr_b,
#                 report_overlap=True,
#             )
#             .assign("frac", lambda df: df.eval("Overlap / (End - Start)"))
#             .subset(lambda df: df["frac"] >= self.fraction)
#             .as_df()
#         )

#         dtype = pd.CategoricalDtype(pr_b.df["Name"].unique())

#         intersection_data = (
#             intersection.set_index("Name")["Name_b"].astype(dtype).rename(self.name)
#         )

#         return intersection_data

#     def _count_intersection(self, pr_b: pr.PyRanges):

#         intersection_data = (
#             self.pr_a.coverage(pr_b)
#             .df.query(f"NumberOverlaps > 0 and FractionOverlaps >= {self.fraction}")
#             .set_index("Name")["NumberOverlaps"]
#             .rename(self.name)
#         )

#         return intersection_data

#     def intersection(self):
#         """
#         Intersect two bed files and return the intersection as a pandas series.

#         Returns:
#             pd.Series: A pandas series containing the intersection.

#         Raises:
#             OSError: Raised if the bed file cannot be read.
#             IndexError: Raised if the bed file is empty.
#             FileNotFoundError: Raised if the bed file cannot be found.
#             StopIteration: Raised if the bed file is empty.
#             AssertionError: Raised if the bed file is empty.

#         """

#         try:

#             pr_b = convert_bed_to_pr(self.b)

#             if self.action == "get":
#                 _intersection = self._get_intersection(pr_b)
#             elif self.action == "count":
#                 _intersection = self._count_intersection(pr_b)

#         except (OSError, IndexError, FileNotFoundError, StopIteration, AssertionError):

#             self.logger.warning(
#                 f"Could not intersect {self.b} using {self.action} method."
#             )
#             _intersection = pd.Series(
#                 data=pd.NA,
#                 index=self.pr_a.df["Name"],
#                 name=self.name,
#                 dtype=object,
#             )

#         return _intersection

#     def __repr__(self):
#         return f"{self.name} intersection"
