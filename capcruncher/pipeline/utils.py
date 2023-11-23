import os
import pathlib
import re
from typing import Dict, List, Union, Literal
from collections import defaultdict
import json
import itertools
import pandas as pd
import pyranges as pr

from capcruncher import utils

from loguru import logger
import snakemake
from snakemake.io import expand, glob_wildcards


def is_on(param: str) -> bool:
    """
    Returns True if parameter in "on" values
    On values:
        - true
        - t
        - on
        - yes
        - y
        - 1
    """
    values = ["true", "t", "on", "yes", "y", "1"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_off(param: str):
    """Returns True if parameter in "off" values"""
    values = ["", "None", "none", "F", "f", "no"]
    if str(param).lower() in values:
        return True
    else:
        return False


def is_none(param: str) -> bool:
    """Returns True if parameter is none"""
    values = ["", "none"]
    if str(param).lower() in values:
        return True
    else:
        return False


def convert_empty_yaml_entry_to_string(param: str) -> str:
    """
    Converts empty yaml entries to string
    """
    if is_none(param):
        return ""
    else:
        return param


def format_config_dict(config: Dict) -> Dict:
    """
    Formats the config dictionary to ensure that all entries are strings.

    """
    for key, value in config.items():
        config[key] = convert_empty_yaml_entry_to_string(value)

    return config


def get_design_matrix(fastqs: List[Union[str, pathlib.Path]]):
    df = pd.DataFrame(fastqs, columns=["fn"])
    df["filename"] = df["fn"].apply(str).str.split(".fastq").str[0]
    df["sample"] = df["filename"].str.extract(r".*/(.*?)_R?[12].fastq.*")
    df["condition"] = df["sample"].str.split(".fastq").str[0].str.split("_").str[-1]

    if df["condition"].isna().any():
        logger.warn(
            "Failed to identify conditions from fastq files. Please format as sample_CONDITION_READ.fastq(.gz)"
        )
        df["condition"].fillna("UNKNOWN")

    return df[["sample_name", "condition"]].drop_duplicates()


def get_bin_sizes(config):
    binsizes = config["analysis"].get("bin_sizes")

    if binsizes is None:
        binsizes = []
    elif isinstance(binsizes, int):
        binsizes = [binsizes]
    elif isinstance(binsizes, str):
        binsizes = [
            int(bs)
            for bs in re.split(
                r"[,;]\s*|\s+", str(config["analysis"].get("bin_sizes", "0"))
            )
        ]
    elif isinstance(binsizes, list):
        binsizes = [int(bs) for bs in binsizes]
    else:
        raise ValueError(
            "Invalid bin size(s). Please specify as int, comma separated string or list of ints."
        )

    return binsizes


def get_blacklist(config):
    if utils.is_valid_bed(config["analysis_optional"].get("blacklist"), verbose=False):
        blacklist = config["analysis_optional"]["blacklist"]
        return blacklist


def has_high_viewpoint_number(viewpoints: str, config: Dict):
    n_viewpoints = pr.read_bed(viewpoints).df.shape[0]
    if n_viewpoints > 500:
        if not config["analysis_optional"].get("force_bigwig_generation", False):
            return True


def can_perform_plotting(config):
    try:
        pass
    except ImportError:
        logger.warning(
            "Plotting capabilities not installed. For plotting please run: pip install capcruncher[plotting]"
        )
        return False

    return utils.is_valid_bed(config["plot"].get("coordinates"), verbose=False)


def can_perform_binning(config):
    perform_binning = False

    try:
        bin_sizes = config["analysis"].get("bin_sizes", None)
        if isinstance(bin_sizes, int) and bin_sizes > 0:
            perform_binning = True
        elif isinstance(bin_sizes, list) and all([int(bs) > 0 for bs in bin_sizes]):
            perform_binning = True

        return perform_binning

    except Exception as e:
        logger.error(e)
        return perform_binning


def group_files_by_regex(files: List, regex: str):
    df = pd.DataFrame(files, columns=["fn"])
    extracted_substrings = df["fn"].astype(str).str.extract(regex)
    df = df.join(extracted_substrings)
    return (
        df.groupby(extracted_substrings.columns.to_list())
        .agg(list)["fn"]
        .rename("files_grouped")
    )


class FastqSamples:
    def __init__(self, design):
        # Expected columns: sample, fq1, fq2
        self._design = design
        self._design = self._design.assign(
            paired=(~self._design[["fq1", "fq2"]].isna().any(axis=1))
        )

    @classmethod
    def from_files(cls, files: List[Union[pathlib.Path, str]]) -> "FastqSamples":
        if not len(files) > 0:
            logger.error("No fastq files found.")
            raise ValueError("No fastq files found.")

        df = pd.DataFrame(files, columns=["fn"])

        df[["sample", "read"]] = (
            df["fn"].apply(str).str.extract("(?!.*/)?(.*).*_R?([12]).fastq.gz")
        )

        df["sample"] = df["sample"].apply(
            lambda p: pathlib.Path(p).name
            if isinstance(p, pathlib.Path)
            else os.path.basename(p)
        )
        df["read"] = "fq" + df["read"]

        df = (
            df.pivot(columns="read", index=["sample"])
            .droplevel(level=0, axis=1)
            .reset_index()
        )

        # Format to check for
        # CONDITION-A_REPLICATE-IDENTIFIER_READNUMBER
        try:
            df[["condition", "replicate"]] = df["sample"].str.split("_", expand=True)
        except ValueError:
            logger.warning("Failed to identify conditions from fastq files.")
            df["condition"] = "UNKNOWN"

        return cls(design=df)

    @property
    def fastq_files(self):
        return sorted([*self._design["fq1"], *self._design["fq2"]])

    @property
    def sample_names_all(self):
        return self._design["sample"].to_list()

    @property
    def translation(self):
        fq_translation = {}
        for sample in self._design.itertuples():
            for read in [1, 2]:
                fq_translation[f"{sample.sample}_{read}.fastq.gz"] = os.path.realpath(
                    str(getattr(sample, f"fq{read}"))
                )

        return fq_translation

    @property
    def design(self):
        return self._design


def validate_blacklist(blacklist):
    """Validate blacklist file."""

    blacklist_ok = True

    if blacklist is None:
        blacklist_ok = False
    elif not os.path.exists(blacklist):
        blacklist_ok = False

    return blacklist_ok


def configure_annotation_parameters(workflow: snakemake.Workflow, config: Dict) -> Dict:
    """Load defaults from annotation_defaults.json and overwrite with the current files"""

    path = pathlib.Path(__file__).absolute()
    defaults_file = path.parent / "workflow" / "data" / "annotation_defaults.json"
    parameters = json.load(open(defaults_file))

    # Overwrite defaults with current options
    parameters["viewpoints"]["fn"] = config["analysis"]["viewpoints"]
    parameters["viewpoints"]["fraction"] = config["analysis_optional"].get(
        "minimum_viewpoint_overlap", parameters["viewpoints"]["fraction"]
    )
    parameters["viewpoints_count"]["fn"] = config["analysis"]["viewpoints"]
    parameters["viewpoints_count"]["fraction"] = config["analysis_optional"].get(
        "minimum_viewpoint_overlap", parameters["viewpoints"]["fraction"]
    )

    # Check if blacklist is valid
    if validate_blacklist(config["analysis"].get("blacklist")):
        parameters["blacklist"] = config["analysis"]["blacklist"]
    else:
        del parameters["blacklist"]

    return parameters


def format_annotation_parameters(*args, **kwargs):
    """Format annotation parameters for use in the shell script."""

    parameters = configure_annotation_parameters(*args)

    flags = {
        "name": "-n",
        "fn": "-b",
        "action": "-a",
        "fraction": "-f",
        "dtype": "-t",
    }

    annotation_args = []
    for annotation, options in parameters.items():
        for option, value in options.items():
            if value is not None:
                annotation_args.append(f"{flags[option]} {value}")

    return " ".join(annotation_args)


def format_priority_chromosome_list(config: Dict):
    """Format priority chromosome list for use in the shell script."""

    priority_chroms = config["analysis_optional"].get("priority_chromosomes", "")

    if not priority_chroms or priority_chroms == "None":
        chromosomes = None
    elif "," in priority_chroms:
        chromosomes = priority_chroms
    elif "viewpoints" in priority_chroms:
        pr_viewpoints = pr.read_bed(config["analysis"]["viewpoints"])
        chromosomes = ",".join(pr_viewpoints.Chromosome)

    return f"--priority-chroms {chromosomes}" if chromosomes else ""


def get_threads_for_annotation(annotation_files_and_params):
    return annotation_files_and_params.count("-b")


def identify_columns_based_on_condition(design: pd.DataFrame):
    condition_args = []

    for group_name, columns in design.groupby("condition").groups.items():
        condition_args.append(f"-c {','.join(str(c) for c in columns)}")

    condition_args_str = " ".join(condition_args)

    return condition_args_str


def validate_custom_filtering(config: Dict):
    custom_filter_stages = config["analysis"].get("custom_filtering", "")
    if not custom_filter_stages:
        cf = ""
    elif not os.path.exists(custom_filter_stages):
        cf = ""
    else:
        cf = f"--custom-filtering {custom_filter_stages}"

    return cf


def get_count_files(wc, perform_binning: bool = False):
    counts = []
    counts.append(
        f"capcruncher_output/interim/pileups/counts_by_restriction_fragment/{wc.sample}.hdf5"
    )

    if perform_binning:
        counts.append(
            f"capcruncher_output/interim/pileups/counts_by_genomic_bin/{wc.sample}.hdf5"
        )

    return counts


def get_normalisation_from_config(wc, config: Dict):
    regions = config["normalisation"]["regions"]

    if regions is not None or isinstance(regions, str):
        if os.path.exists(regions):  # noqa: E714
            return f"--normalisation region --normalisation-regions {regions}"
    return "--normalisation n_cis"


def get_fastq_basename(wildcards, fastq_samples: FastqSamples, **kwargs):
    return pathlib.Path(
        fastq_samples.translation[f"{wildcards.sample}_{wildcards.read}.fastq.gz"]
    ).stem.replace(".fastq", "")


def get_files_to_plot(
    wc,
    design: pd.DataFrame,
    assay: Literal["capture", "tri", "tiled"],
    sample_names: List[str],
    summary_methods: List[str],
    compare_samples: bool = False,
):
    files = {
        "bigwigs": [],
        "subtractions": [],
        "bigwigs_collection": [],
        "heatmaps": [],
    }

    if assay == "tiled":
        files["heatmaps"].extend(
            expand(
                "capcruncher_output/results/{sample}/{sample}.hdf5",
                sample=sample_names,
            )
        )
        return files

    if compare_samples:
        bigwigs_comparison = expand(
            "capcruncher_output/results/comparisons/bigwigs/{comparison}.{method}-subtraction.{{viewpoint}}.bigWig",
            comparison=[
                f"{a}-{b}"
                for a, b in itertools.permutations(design["condition"].unique(), 2)
            ],
            method=summary_methods,
        )

        files["subtractions"].extend(bigwigs_comparison)

    bigwigs = expand(
        "capcruncher_output/results/{sample}/bigwigs/norm/{sample}_{{viewpoint}}.bigWig",
        sample=sample_names,
    )

    files["bigwigs"].extend(bigwigs)

    return files


def get_plotting_coordinates(wc, config: Dict):
    plot_coords = config["plot"].get("coordinates", None)

    if plot_coords and pathlib.Path(plot_coords).exists():
        df = pd.read_table(
            plot_coords, names=["chrom", "start", "end", "name"], header=None
        )
        df = df.query("name.str.contains(@wc.viewpoint)").iloc[0]

    else:
        df = pd.read_table(
            config["analysis"]["viewpoints"],
            names=["chrom", "start", "end", "name"],
            header=None,
        )

        df = df.query("name == @wc.viewpoint").iloc[0]

    return f"{df.chrom}:{df.start}-{df.end}"


def get_pileups(
    assay: Literal["capture", "tri", "tiled"],
    design: pd.DataFrame,
    samples_aggregate: bool,
    samples_compare: bool,
    sample_names: List[str],
    summary_methods: List[str],
    viewpoints: List[str],
) -> list[str]:
    bigwigs = []
    if assay in ["capture", "tri"]:
        bigwigs.extend(
            expand(
                "capcruncher_output/results/{sample}/bigwigs/{norm}/{sample}_{viewpoint}.bigWig",
                sample=sample_names,
                norm=["raw", "norm"],
                viewpoint=viewpoints,
            ),
        )

        if samples_aggregate:
            bigwigs.extend(
                expand(
                    "capcruncher_output/results/comparisons/bigwigs/{group}.{method}-summary.{viewpoint}.bigWig",
                    group=design["condition"].unique(),
                    method=summary_methods,
                    viewpoint=viewpoints,
                ),
            )

        if samples_compare:
            bigwigs.extend(
                expand(
                    "capcruncher_output/results/comparisons/bigwigs/{comparison}.{method}-subtraction.{viewpoint}.bigWig",
                    comparison=[
                        f"{a}-{b}"
                        for a, b in itertools.permutations(
                            design["condition"].unique(), 2
                        )
                    ],
                    method=summary_methods,
                    viewpoint=viewpoints,
                ),
            )

    elif assay == "tiled":
        pass

    return bigwigs
