import os
import sys
from click.types import STRING
import pandas as pd
import numpy as np
from pybedtools import BedTool
import cooler
import h5py
from joblib import Parallel, delayed
from typing import Union, Literal
from natsort import natsorted, natsort_key
from capcruncher.utils import split_intervals_on_chrom, intersect_bins
import itertools


def get_capture_coords(viewpoint_file: str, viewpoint_name: str):
    df_viewpoints = BedTool(viewpoint_file).to_dataframe()
    df_viewpoints = df_viewpoints.query(f'name == "{viewpoint_name}"')

    try:
        viewpoints = [row for index, row in df_viewpoints.iterrows()]
    except IndexError:
        print("Oligo name cannot be found within viewpoints")
        viewpoints = None

    return viewpoints


def get_capture_bins(bins, viewpoint_chrom, viewpoint_start, viewpoint_end):

    return [
        int(b)
        for b in bins.query(
            f'chrom == "{viewpoint_chrom}" and start >= {viewpoint_start} and end <= {viewpoint_end}'
        )["name"]
    ]


def create_cooler_cc(
    output_prefix: str,
    bins: pd.DataFrame,
    pixels: pd.DataFrame,
    capture_name: str,
    capture_viewpoints: os.PathLike,
    capture_bins: Union[int, list] = None,
    suffix=None,
    **cooler_kwargs,
) -> os.PathLike:
    """
    Creates a cooler hdf5 file or cooler formatted group within a hdf5 file.

    Args:
     output_prefix (str): Output path for hdf5 file. If this already exists, will append a new group to the file.
     bins (pd.DataFrame): DataFrame containing the genomic coordinates of all bins in the pixels table.
     pixels (pd.DataFrame): DataFrame with columns: bin1_id, bin2_id, count.
     capture_name (str): Name of capture probe to store.
     capture_viewpoints (os.PathLike): Path to capture viewpoints used for the analysis.
     capture_bins (Union[int, list], optional): Bins containing capture viewpoints. Can be determined from viewpoints if not supplied. Defaults to None.
     suffix (str, optional): Suffix to append before the .hdf5 file extension. Defaults to None.

    Raises:
     ValueError: Capture name must exactly match the name of a supplied capture viewpoint.

    Returns:
     os.PathLike: Path of cooler hdf5 file.
    """

    # Gets capture coordinates
    capture_coords = get_capture_coords(capture_viewpoints, capture_name)

    # Make sure capture coordinates are returned correctly, if not, error.
    if capture_coords is None:
        raise ValueError(f"Incorrect capture name specified: {capture_name}.")

    # If capture bins not provided get them using the coordinates.
    if not capture_bins:
        capture_bins = list(
            itertools.chain.from_iterable(
                [
                    get_capture_bins(bins, c["chrom"], c["start"], c["end"])
                    for c in capture_coords
                ]
            )
        )

    # Need to store bins as a list so make sure its not just a single int.
    elif isinstance(capture_bins, int):
        capture_bins = [
            int(capture_bins),
        ]

    # The cooler.create_cooler function will not accept np.arrays so must convert to python list
    elif isinstance(capture_bins, (np.array, pd.Series)):
        capture_bins = [int(x) for x in capture_bins]

    # Get the number of cis interactions, required for normalisation.
    bins_cis = bins.loc[
        lambda df: df["chrom"].isin([c["chrom"] for c in capture_coords])
    ]["name"].loc[lambda ser: ~ser.isin(capture_bins)]

    pixels_cis = pixels.loc[
        lambda df: (df["bin1_id"].isin(bins_cis)) | (df["bin2_id"].isin(bins_cis))
    ]
    n_cis_interactions = pixels_cis["count"].sum()

    # Metadata for cooler file.
    metadata = {
        "capture_bins": capture_bins,
        "capture_name": capture_name,
        "capture_chrom": [c["chrom"] for c in capture_coords],
        "capture_coords": [
            f'{c["chrom"]}:{c["start"]}-{c["end"]}' for c in capture_coords
        ],
        "n_cis_interactions": int(n_cis_interactions),
    }

    if os.path.exists(
        output_prefix
    ):  # Will append to a prexisting file if one is supplied
        append_to_file = True
        cooler_fn = f"{output_prefix}::/{capture_name}"
    else:
        append_to_file = False
        cooler_fn = f"{output_prefix.replace('.hdf5', '')}.{capture_name}{'.' + suffix if suffix else ''}.hdf5"

    cooler.create_cooler(
        cooler_fn,
        bins=bins,
        pixels=pixels,
        metadata=metadata,
        mode="w" if not append_to_file else "a",
        **cooler_kwargs,
    )

    return cooler_fn


class GenomicBinner:
    """
    Provides a conversion table for converting two sets of bins.

    Attributes:
     chromsizes (pd.Series): Series indexed by chromosome name containg chromosome sizes in bp
     fragments (pd.DataFrame): DataFrame containing bins to convert to equal genomic intervals
     binsize (int): Genomic bin size
     min_overlap (float): Minimum degree of intersection to define an overlap.
     n_cores (int): Number of cores to use for bin intersection.

    """

    def __init__(
        self,
        chromsizes: Union[os.PathLike, pd.DataFrame, pd.Series],
        fragments: pd.DataFrame,
        binsize: int = 5000,
        n_cores: int = 8,
        method: Literal["midpoint", "overlap"] = "midpoint",
        min_overlap: float = 0.2,
    ):
        """
        Args:
         chromsizes (Union[os.PathLike, pd.DataFrame, pd.Series]): Series indexed by chromosome name containg chromosome sizes in bp
         fragments (pd.DataFrame): DataFrame containing bins to convert to equal genomic intervals
         binsize (int, optional): Genomic window size. Defaults to 5000.
         n_cores (int, optional): Number of cores to use for bin intersection.. Defaults to 8.
         min_overlap (float, optional): Minimum degree of intersection to define an overlap.Only used for "overlap" method. Defaults to 0.2.
        """

        if not method in ["midpoint", "overlap"]:
            raise ValueError('Method should be either "midpoint" or "overlap"')

        self.method = method
        self.chromsizes = self._format_chromsizes(chromsizes)
        self.fragments = self._format_fragments(fragments)
        self.binsize = binsize
        self.min_overlap = min_overlap

        self.bins_genomic = self._get_bins()
        self._bin_conversion_table = None
        self.n_cores = n_cores

    def _get_bins(self):
        return (
            cooler.util.make_bintable(chromsizes=self.chromsizes, binsize=self.binsize)
            .reset_index()
            .rename(columns={"index": "name"})[["chrom", "start", "end", "name"]]
            .assign(
                start=lambda df: df["start"].astype(int),
                end=lambda df: df["end"].astype(int),
            )
        )

    def _format_chromsizes(self, chromsizes):

        _chromsizes = pd.Series(dtype=np.int64)
        if isinstance(chromsizes, str):
            _chromsizes = pd.read_csv(
                chromsizes,
                sep="\t",
                header=None,
                names=["chrom", "size"],
                usecols=[0, 1],
                index_col=0,
            )["size"].sort_index(key=natsort_key)

        elif isinstance(chromsizes, pd.DataFrame):
            if chromsizes.index.astype(str).str.contains("^chr.*"):
                _chromsizes = chromsizes.iloc[:, 0]

        elif isinstance(chromsizes, pd.Series):
            _chromsizes = chromsizes

        if not _chromsizes.empty:
            return _chromsizes
        else:
            raise ValueError("Chromsizes supplied in the wrong format")

    def _natsort_dataframe(self, df, column):

        df_by_key = {k: df for k, df in df.groupby(column)}

        _df = pd.DataFrame()
        for k in natsorted(df_by_key):
            if _df is not None:
                _df = pd.concat([_df, df_by_key[k]])
            else:
                _df = df_by_key[k]

        return _df

    def _format_fragments(self, fragments):

        # Read the fragments
        if isinstance(fragments, str):
            _fragments = pd.read_csv(
                fragments,
                sep="\t",
                index_col=0,
                header=None,
                names=["chrom", "start", "end", "name"],
            )
        elif isinstance(fragments, pd.DataFrame):
            _fragments = fragments

            if not "name" in _fragments.columns:
                _fragments = _fragments.reset_index().rename(columns={"index": "name"})[
                    ["chrom", "start", "end", "name"]
                ]

        # Adjust fragments based on method
        if self.method == "midpoint":
            lengths = _fragments["end"] - fragments["start"]
            midpoint = _fragments["start"] + (lengths // 2)

            _fragments["start"] = midpoint
            _fragments["end"] = midpoint + 1

        return self._natsort_dataframe(_fragments, "chrom")

    def _get_bin_conversion_table(self) -> pd.DataFrame:

        # Split bins by chromosome
        bins_genomic_by_chrom = split_intervals_on_chrom(self.bins_genomic)
        bins_fragments_by_chrom = split_intervals_on_chrom(self.fragments)

        # Identify shared chromosomes (all should be represented)
        shared_chroms = set(bins_fragments_by_chrom) & set(bins_genomic_by_chrom)

        # Perform intersection of bins by chromosome
        bins_intersections = Parallel(n_jobs=self.n_cores)(
            delayed(intersect_bins)(
                bins_fragments_by_chrom[chrom],
                bins_genomic_by_chrom[chrom],
                loj=True,
                sorted=True,
                f=self.min_overlap if self.method == "overlap" else 1e-9,
            )
            for chrom in natsorted(shared_chroms)
        )

        # Concatenate intersected dataframes and replace labels for clarity
        df_bins_intersections = pd.concat(bins_intersections, ignore_index=True).rename(
            columns=lambda c: c.replace("_1", "_fragment").replace("_2", "_bin")
        )

        # # Calculate overlap fraction
        # df_bins_intersections["overlap_fraction"] = df_bins_intersections["overlap"] / (
        #     df_bins_intersections["end_fragment"]
        #     - df_bins_intersections["start_fragment"]
        # )

        return df_bins_intersections

    @property
    def bins(self) -> pd.DataFrame:
        """
        Equal genomic bins.

        Returns:
         pd.DataFrame: DataFrame in bed format.
        """
        return self.bins_genomic

    @property
    def bin_conversion_table(self) -> pd.DataFrame:
        """

        Returns:
         pd.DataFrame: Conversion table containing coordinates and ids of intersecting bins.
        """

        if self._bin_conversion_table is not None:
            return self._bin_conversion_table
        else:
            self._bin_conversion_table = self._get_bin_conversion_table()
            return self._bin_conversion_table


class CoolerBinner:
    """Bins a cooler file into equal genomic intervals.

    Attributes:
     cooler: (cooler.Cooler): Cooler instance to bin.
     binner (capcruncher.storeage.GenomicBinner): Binner class to generate bin conversion tables
     binsize (int): Genomic bin size
     scale_factor (int): Scaling factor for normalising interaction counts.
     n_cis_interactions (int): Number of cis interactions with the capture bins.
     n_cores (int): Number of cores to use for binning.

    """

    def __init__(
        self,
        cooler_fn: os.PathLike,
        binsize: int = None,
        n_cores: int = 8,
        binner: GenomicBinner = None,
    ):
        """
        Args:
         cooler_fn (os.PathLike): Path to cooler to bin. A cooler group can be specified with FN_PATH.hdf5::/PATH_TO_GROUP.
         binsize (int, optional): Genomic binsize. Defaults to None.
         n_cores (int, optional): Number of cores to use for binning. Defaults to 8.
         binner (capcruncher.storage.GenomicBinner, optional): Binner object to produce conversion tables.
                                                                 Can be initialised and provided so that binning is not repeated.
                                                                 Defaults to None.
        """

        self.cooler = cooler.Cooler(cooler_fn)
        self.bins_fragments = self.cooler.bins()[:]

        self.binner = binner or GenomicBinner(
            chromsizes=self.cooler.chromsizes,
            fragments=self.bins_fragments,
            binsize=binsize,
        )

        self.binsize = self.binner.binsize
        self.n_cores = n_cores

        self._bin_conversion_table = None
        self._pixel_conversion_table = None
        self._pixels = None

        self.n_cis_interactions = self.cooler.info["metadata"]["n_cis_interactions"]

    def _get_pixel_conversion_table(self):

        pixels = self.cooler.pixels()[:]
        pixels_conv = (
            pixels.merge(
                self.bin_conversion_table[["name_bin", "name_fragment"]].add_suffix(
                    "_1"
                ),
                left_on="bin1_id",
                right_on="name_fragment_1",
            )
            .merge(
                self.bin_conversion_table[["name_bin", "name_fragment"]].add_suffix(
                    "_2"
                ),
                left_on="bin2_id",
                right_on="name_fragment_2",
            )
            .drop(columns=["name_fragment_1", "name_fragment_2"])
        )

        return pixels_conv

    def _get_pixels(self):

        df_pixels = (
            self.pixel_conversion_table.groupby(["name_bin_1", "name_bin_2"])["count"]
            .sum()
            .to_frame()
            .reset_index()
        )

        df_pixels.columns = ["bin1_id", "bin2_id", "count"]

        # Swap bins over if not correct
        df_pixels["bin1_id_corrected"] = np.where(
            df_pixels["bin1_id"] > df_pixels["bin2_id"],
            df_pixels["bin2_id"],
            df_pixels["bin1_id"],
        )
        df_pixels["bin2_id_corrected"] = np.where(
            df_pixels["bin1_id"] > df_pixels["bin2_id"],
            df_pixels["bin1_id"],
            df_pixels["bin2_id"],
        )
        # df_pixels = df_pixels.loc[lambda df: df["bin1_id"] != df["bin2_id"]]
        df_pixels = df_pixels.loc[
            :, ["bin1_id_corrected", "bin2_id_corrected", "count"]
        ].rename(columns=lambda col: col.replace("_corrected", ""))

        return df_pixels

    @property
    def bins(self) -> pd.DataFrame:
        """
        Returns:
         pd.DataFrame: Even genomic bins of a specified binsize.
        """
        return self.binner.bins

    @property
    def bin_conversion_table(self) -> pd.DataFrame:
        """
        Returns:
         pd.DataFrame: Conversion table containing coordinates and ids of intersecting bins.
        """
        return self.binner.bin_conversion_table

    @property
    def capture_bins(self):
        """
        Returns:
         pd.DataFrame: Capture bins converted to the new even genomic bin format.
        """
        capture_frags = self.cooler.info["metadata"]["capture_bins"]
        return self.bin_conversion_table.loc[
            lambda df: df["name_fragment"].isin(capture_frags)
        ]["name_bin"].values

    @property
    def pixel_conversion_table(self):
        """
        Returns:
         pd.DataFrame: Conversion table to convert old binning scheme to the new.
        """
        if self._pixel_conversion_table is not None:
            return self._pixel_conversion_table
        else:
            self._pixel_conversion_table = self._get_pixel_conversion_table()
            return self._pixel_conversion_table

    @property
    def pixels(self):
        """
        Returns:
         pd.DataFrame: Pixels (interaction counts) converted to the new binning scheme.
        """
        if self._pixels is not None:
            return self._pixels
        else:
            self._pixels = self._get_pixels()
            return self._pixels

    def normalise_pixels(
        self,
        n_fragment_correction: bool = True,
        n_interaction_correction: bool = True,
        scale_factor: int = 1e6,
    ):
        """Normalises pixels (interactions).

        Normalises pixels according to the number of restriction fragments per bin and the number of cis interactions.
        If both normalisation options are selected, will also provide a dual normalised column.

        Args:
         n_fragment_correction (bool, optional): Updates the pixels DataFrame with counts corrected for the number of
                                                 restriction fragments per bin. Defaults to True.
         n_interaction_correction (bool, optional): Updates the pixels DataFrame with counts corrected for the number of
                                                    cis interactions. Defaults to True.
         scale_factor (int, optional): Scaling factor for n_interaction_correction. Defaults to 1e6.
        """

        if n_fragment_correction:
            df_nrf = (
                self.bin_conversion_table.groupby("name_bin").size().to_frame("n_rf")
            )
            pixels = self.pixels.merge(
                df_nrf.add_prefix("bin1_"), left_on="bin1_id", right_index=True
            ).merge(df_nrf.add_prefix("bin2_"), left_on="bin2_id", right_index=True)

            self.pixels["count_n_rf_norm"] = self.pixels["count"] / (
                pixels["bin1_n_rf"] * pixels["bin2_n_rf"]
            )

        if n_interaction_correction:
            self.pixels["count_n_interactions_norm"] = (
                self.pixels["count"] / self.n_cis_interactions
            ) * scale_factor

        if n_fragment_correction and n_interaction_correction:

            self.pixels["count_n_rf_n_interactions_norm"] = (
                self.pixels["count_n_rf_norm"] / self.n_cis_interactions
            ) * scale_factor

    def to_cooler(self, store, normalise=False, **normalise_options):

        capture_bins = self.capture_bins
        capture_name = self.cooler.info["metadata"]["capture_name"]
        capture_coords = self.cooler.info["metadata"]["capture_coords"]
        capture_chrom = self.cooler.info["metadata"]["capture_chrom"]

        metadata = {
            "capture_bins": [int(x) for x in self.capture_bins],
            "capture_name": capture_name,
            "capture_coords": capture_coords,
            "capture_chrom": capture_chrom,
            "n_cis_interactions": self.n_cis_interactions,
        }

        if normalise:
            self.normalise_pixels(**normalise_options)

        if os.path.exists(store):  # Will append to a prexisting file if one is supplied
            cooler_fn = f"{store}::/{capture_name}/resolutions/{self.binsize}"
        else:
            cooler_fn = (
                f"{store.replace('.hdf5', '')}.{capture_name}.{self.binsize}.hdf5"
            )

        cooler.create_cooler(
            cooler_fn,
            bins=self.bins,
            pixels=self.pixels,
            metadata=metadata,
            mode="w" if not os.path.exists(store) else "a",
            columns=self.pixels.columns[2:],
        )

        return cooler_fn


def link_bins(clr: os.PathLike):
    """Reduces cooler storage space by linking "bins" table.

     All of the cooler "bins" tables containing the genomic coordinates of each bin
     are identical for all cooler files of the same resoultion. As cooler.create_cooler
     generates a new bins table for each cooler, this leads to a high degree of duplication.

     This function hard links the bins tables for a given resolution to reduce the degree of duplication.

    Args:
     clr (os.PathLike): Path to cooler hdf5 produced by the merge command.
    """

    with h5py.File(clr, "a") as f:

        # Get all captures stored
        captures = list(f.keys())

        # Get all resolutions stored
        resolutions_group = f[captures[0]].get("resolutions")
        resolutions = list(resolutions_group.keys()) if resolutions_group else None

        for capture in captures[1:]:

            # Delete currenly stored bins group and replace with link to first capture "bins" group
            del f[capture]["bins"]
            f[capture]["bins"] = f[captures[0]]["bins"]

            if resolutions:
                for resolution in resolutions:
                    # Repeat for resolutions i.e. binned coolers
                    del f[capture]["resolutions"][resolution]["bins"]
                    f[capture]["resolutions"][resolution]["bins"] = f[captures[0]][
                        "resolutions"
                    ][resolution]["bins"]
