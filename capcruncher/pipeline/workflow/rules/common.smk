import math
import os
import pathlib
import re
import shutil
import sys
from typing import Literal
from typing import List


def scale_resource(base_value: int, attempt: int = 1) -> int:
    return max(1, math.ceil(base_value * 2 ** (int(attempt) - 1) * SCALE_RESOURCES))


def scale_memory(base_gb: int, attempt: int = 1) -> str:
    return f"{scale_resource(base_gb, attempt)}G"


def scale_thread_memory(base_gb: int, threads: int, attempt: int = 1) -> str:
    return scale_memory(base_gb * threads, attempt)


def default_fastq_split_method() -> str:
    configured_method = config.get("split", {}).get("method")
    if configured_method:
        return configured_method

    if sys.platform == "darwin" and shutil.which("gsplit") is None:
        return "python"

    return "unix"


def copy_workflow_log(log_value, destination):
    if isinstance(log_value, (str, bytes, os.PathLike)):
        shutil.copyfile(log_value, destination)
        return

    for log_path in reversed(list(log_value)):
        if pathlib.Path(log_path).exists():
            shutil.copyfile(log_path, destination)
            return


def get_split_1_parts(wildcards):
    outdir = checkpoints.split.get(**wildcards).output[0]
    fq_files = pathlib.Path(outdir).glob("*.fastq.gz")
    return sorted(
        {
            int(re.search(r"part(\d+)", f.name).group(1))
            for f in fq_files
            if re.search(r"part(\d+)", f.name)
        }
    )


def get_pickles(wc):
    return expand(
        "capcruncher_output/interim/fastq/deduplicated/{{sample}}/{{sample}}_{part}.pkl",
        part=get_split_1_parts(wc),
    )


def get_fastq_split_1(wildcards):
    return {
        f"fq{read}": expand(
            "capcruncher_output/interim/fastq/split/{{sample}}/{{sample}}_part{part}_{read}.fastq.gz",
            part=get_split_1_parts(wildcards),
            read=[read],
        )
        for read in ["1", "2"]
    }


def get_deduplicated_fastq_pair(wildcards):
    input_dir = checkpoints.deduplication.get(**wildcards).output[0]
    fq = {
        f"fq{read}": f"{input_dir.rstrip('/')}/{wildcards.sample}_part{wildcards.part}_{read}.fastq.gz"
        for read in ["1", "2"]
    }

    if pathlib.Path(fq["fq1"]).exists() and pathlib.Path(fq["fq2"]).exists():
        return fq
    return {"fq1": [], "fq2": []}


def get_flashed_fastq(wildcards):
    return [
        f"capcruncher_output/interim/fastq/flashed/{wildcards.sample}/{wildcards.sample}_part{part}.extendedFrags.fastq.gz"
        for part in get_split_1_parts(wildcards)
    ]


def get_pe_fastq(wildcards):
    return [
        f"capcruncher_output/interim/fastq/flashed/{wildcards.sample}/{wildcards.sample}_part{part}.notCombined_{read}.fastq.gz"
        for part in get_split_1_parts(wildcards)
        for read in ["1", "2"]
    ]


def get_rebalanced_parts(
    wildcards, combined: Literal["flashed", "pe"] = None, **kwargs
):
    combined = combined or wildcards.combined
    parts = {}
    outdirs = {
        "flashed": checkpoints.rebalance_partitions_combined.get(
            **{**wildcards, **kwargs}
        ).output[0],
        "pe": checkpoints.rebalance_partitions_pe.get(**{**wildcards, **kwargs}).output[
            0
        ],
    }

    for combined_type in ["flashed", "pe"]:
        fq_files = pathlib.Path(outdirs[combined_type]).glob("*.fastq.gz")
        parts[combined_type] = sorted(
            {
                int(re.search(r"part(\d+)", f.name).group(1))
                for f in fq_files
                if re.search(r"part(\d+)", f.name)
            }
        )

    if combined == "flashed":
        return parts["flashed"]
    return parts["pe"]


def get_rebalanced_fastq_combined(wc):
    checkpoint_output = checkpoints.rebalance_partitions_combined.get(**wc).output[0]
    return f"capcruncher_output/interim/fastq/rebalanced/{wc.sample}/flashed/{wc.sample}_part{wc.part}_flashed_1.fastq.gz"


def get_rebalanced_fastq_pe(wc):
    checkpoint_output = checkpoints.rebalance_partitions_pe.get(**wc).output[0]
    return {
        "pe1": f"capcruncher_output/interim/fastq/rebalanced/{wc.sample}/pe/{wc.sample}_part{wc.part}_pe_1.fastq.gz",
        "pe2": f"capcruncher_output/interim/fastq/rebalanced/{wc.sample}/pe/{wc.sample}_part{wc.part}_pe_2.fastq.gz",
    }


def get_deduplicated_fastq(wc):
    checkpoint_output = checkpoints.deduplication.get(sample=wc.sample).output[0]
    return {
        "fq1": f"capcruncher_output/interim/fastq/deduplicated/{wc.sample}/{wc.sample}_part{wc.part}_1.fastq.gz",
        "fq2": f"capcruncher_output/interim/fastq/deduplicated/{wc.sample}/{wc.sample}_part{wc.part}_2.fastq.gz",
    }


def separate_pe_fastq(wc):
    return {
        1: expand(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_{read}.fastq.gz",
            sample=wc.sample,
            part=get_split_1_parts(wc),
            read=["1"],
        ),
        2: expand(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_{read}.fastq.gz",
            sample=wc.sample,
            part=get_split_1_parts(wc),
            read=["2"],
        ),
    }


def get_rebalanced_bam(wildcards):
    bam = []
    for combined_type in ["flashed", "pe"]:
        for part in get_rebalanced_parts(wildcards, combined_type):
            bam.append(
                f"capcruncher_output/interim/aligned/{wildcards.sample}/{wildcards.sample}_part{part}_{combined_type}.sorted.bam"
            )

    return bam


def get_filtered_slices(wildcards):
    slices = {}
    for combined_type in ["flashed", "pe"]:
        parts = get_rebalanced_parts(wildcards, combined=combined_type)
        slices[combined_type] = [
            f"capcruncher_output/interim/filtering/initial/{wildcards.sample}/{wildcards.sample}_part{part}_{combined_type}.slices.parquet"
            for part in parts
        ]
    return slices


def get_annotated_slices(wildcards):
    slices = {}
    for combined_type in ["flashed", "pe"]:
        parts = get_rebalanced_parts(wildcards, combined=combined_type)
        slices[combined_type] = [
            f"capcruncher_output/interim/annotate/{wildcards.sample}/{wildcards.sample}_part{part}_{combined_type}.parquet"
            for part in parts
        ]
    return [*slices["flashed"], *slices["pe"]]


def get_mem(wildcards, threads, attempt=1):
    return scale_thread_memory(3, threads, attempt)


def get_outdir(wildcards, output):
    return str(pathlib.Path(output[0]).parent)


def get_digestion_statistics(wc, sample_names: List[str]):
    stat_files = []
    for sample in sample_names:
        for combined in ["flashed", "pe"]:
            for part in get_rebalanced_parts(wc, combined=combined, sample=sample):
                stat_files.append(
                    f"capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_{combined}.json"
                )

    return stat_files


def get_filtering_statistics(wc, sample_names: List[str]):
    stat_files = []
    for sample in sample_names:
        for combined in ["flashed", "pe"]:
            for part in get_rebalanced_parts(wc, combined=combined, sample=sample):
                stat_files.append(
                    f"capcruncher_output/interim/statistics/filtering/data/{sample}_part{part}_{combined}.json"
                )

    return stat_files
