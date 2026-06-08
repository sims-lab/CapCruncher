from __future__ import annotations

import itertools
import json
import os
import pathlib
import re
from collections.abc import Sequence
from typing import Self

import pandas as pd
import pyranges1 as pr
import yaml
from loguru import logger
from snakemake.io import expand

from capcruncher import utils
from capcruncher.pipeline.validation import (
    format_pipeline_config,
    is_none,
)
from capcruncher.types import Assay


def convert_empty_yaml_entry_to_string(param: str) -> str:
    """
    Converts empty yaml entries to string
    """
    if is_none(param):
        return ""
    else:
        return param


def _genome_profiles_dir() -> pathlib.Path:
    xdg = os.environ.get("XDG_CONFIG_HOME")
    base = (
        pathlib.Path(xdg).expanduser() if xdg else pathlib.Path.home() / ".capcruncher"
    )
    return base / "genomes"


def load_genome_profile(name: str) -> dict:
    profiles_dir = _genome_profiles_dir()
    profile_path = profiles_dir / f"{name}.yml"
    if not profile_path.exists():
        available = (
            [p.stem for p in sorted(profiles_dir.glob("*.yml"))]
            if profiles_dir.exists()
            else []
        )
        hint = (
            f"Available profiles: {', '.join(available)}"
            if available
            else "No profiles found. Run `capcruncher genome add` to create one."
        )
        raise ValueError(f"Genome profile '{name}' not found. {hint}")
    return yaml.safe_load(profile_path.read_text())


def format_config_dict(config: dict) -> dict:
    """
    Normalise and validate the pipeline config in place.

    Resolves genome profiles before Pydantic validation so that
    ``genome: {profile: hg38}`` expands to the full genome block.
    """
    genome_section = config.get("genome", {})
    if isinstance(genome_section, dict) and "profile" in genome_section:
        profile_name = genome_section.pop("profile")
        profile_data = load_genome_profile(profile_name)
        config["genome"] = {**profile_data, **genome_section}

    return format_pipeline_config(config)


def infer_design_from_fastqs(
    fastqs: Sequence[str | pathlib.Path],
    condition_pattern: str | None = None,
) -> pd.DataFrame:
    """Infer a design matrix from FASTQ filenames.

    Expected convention: <CONDITION>_<REPLICATE>_R[12].fastq[.gz]
    condition_pattern: optional regex with a named group ``condition`` to
    override the default rsplit-based extraction.
    """
    df = pd.DataFrame(fastqs, columns=["fn"])
    df["filename"] = df["fn"].apply(lambda fn: pathlib.Path(fn).name)
    df["sample"] = df["filename"].str.extract(r"(.+)_R?[12]\.fastq(?:\.gz)?$")

    if condition_pattern:
        extracted = df["sample"].str.extract(condition_pattern)
        df["condition"] = (
            extracted["condition"] if "condition" in extracted.columns else pd.NA
        )
    else:
        df["condition"] = df["sample"].str.rsplit("_", n=1).str[0]

    df["replicate"] = df["sample"].str.rsplit("_", n=1).str[1]

    if df["condition"].isna().any():
        logger.warning(
            "Could not infer conditions from FASTQ names. "
            "Expected <CONDITION>_<REPLICATE>_R[12].fastq[.gz]. "
            "Setting condition to UNKNOWN."
        )
        df["condition"] = df["condition"].fillna("UNKNOWN")

    return df[["sample", "condition", "replicate"]].drop_duplicates(subset=["sample"])


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


def has_high_viewpoint_number(viewpoints: str, config: dict) -> bool | None:
    n_viewpoints = pr.read_bed(pathlib.Path(viewpoints)).shape[0]
    if n_viewpoints > 500:
        if not config["analysis_optional"].get("force_bigwig_generation", False):
            return True


def can_perform_plotting(config):
    try:
        import plotnado  # noqa: F401
    except ImportError:
        logger.warning(
            "Plotting capabilities not installed. For plotting please run: pip install capcruncher[plot]"
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


def group_files_by_regex(files: Sequence, regex: str) -> pd.Series:
    df = pd.DataFrame(files, columns=["fn"])
    extracted_substrings = df["fn"].astype(str).str.extract(regex)
    df = df.join(extracted_substrings)
    return (
        df.groupby(extracted_substrings.columns.to_list())
        .agg(list)["fn"]
        .rename("files_grouped")
    )


def is_valid_viewpoint_name(name: str):
    return re.match(r"^[A-Za-z0-9_\-]+$", name)


class FastqSamples:
    def __init__(self, design):
        # Expected columns: sample, fq1, fq2
        self._design = design
        self._design = self._design.assign(
            paired=(~self._design[["fq1", "fq2"]].isna().any(axis=1))
        )

    @classmethod
    def from_files(cls, files: Sequence[pathlib.Path | str]) -> Self:
        if not len(files) > 0:
            logger.error("No fastq files found.")
            raise ValueError("No fastq files found.")

        df = pd.DataFrame(files, columns=["fn"])

        df[["sample", "read"]] = (
            df["fn"].apply(str).str.extract("(?!.*/)?(.*).*_R?([12]).fastq.gz")
        )

        df["sample"] = df["sample"].apply(
            lambda p: (
                pathlib.Path(p).name
                if isinstance(p, pathlib.Path)
                else os.path.basename(p)
            )
        )
        df["read"] = "fq" + df["read"]

        df = (
            df.pivot(columns="read", index=["sample"])
            .droplevel(level=0, axis=1)
            .reset_index()
        )

        design_info = infer_design_from_fastqs(files)
        df = df.merge(
            design_info[["sample", "condition", "replicate"]], on="sample", how="left"
        )

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


def configure_annotation_parameters(workflow, config: dict) -> dict:
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
    for _, options in parameters.items():
        for option, value in options.items():
            if value is not None:
                annotation_args.append(f"{flags[option]} {value}")

    return " ".join(annotation_args)


def format_priority_chromosome_list(config: dict):
    """Format priority chromosome list for use in the shell script."""

    priority_chroms = config["analysis_optional"].get("priority_chromosomes", "")
    chromosomes = None

    if not priority_chroms or priority_chroms == "None":
        chromosomes = None
    elif "," in priority_chroms:
        chromosomes = priority_chroms
    elif "viewpoints" in priority_chroms:
        pr_viewpoints = pr.read_bed(pathlib.Path(config["analysis"]["viewpoints"]))
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


def validate_filter_profile(config: dict):
    filter_profile = config["analysis"].get("filter_profile", "")
    if not filter_profile:
        cf = ""
    elif not os.path.exists(filter_profile):
        cf = ""
    else:
        cf = f"--filter-profile {filter_profile}"

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


def get_normalisation_from_config(wc, config: dict):
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
    assay: Assay,
    sample_names: list[str],
    summary_methods: list[str],
    compare_samples: bool = False,
):
    files = {
        "bigwigs": [],
        "subtractions": [],
        "bigwigs_collection": [],
        "heatmaps": [],
    }

    if assay == Assay.TILED:
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
                f"{a}_vs_{b}"
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


def get_plotting_coordinates(wc, config: dict):
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
    assay: Assay,
    design: pd.DataFrame,
    samples_aggregate: bool,
    samples_compare: bool,
    sample_names: list[str],
    summary_methods: list[str],
    viewpoints: list[str],
) -> list[str]:
    bigwigs = []
    if assay in {Assay.CAPTURE, Assay.TRI}:
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
                        f"{a}_vs_{b}"
                        for a, b in itertools.permutations(
                            design["condition"].unique(), 2
                        )
                    ],
                    method=summary_methods,
                    viewpoint=viewpoints,
                ),
            )

    elif assay == Assay.TILED:
        pass

    return bigwigs
