import glob
import os
import re
import shutil
import subprocess
import sys
from collections.abc import Sequence
from multiprocessing import SimpleQueue
from pathlib import Path
from typing import Any, Literal

from joblib import Parallel, delayed
from loguru import logger

PLATFORM = sys.platform
type FilePath = str | os.PathLike[str]


def run_unix_split(
    fn: FilePath,
    n_reads: int,
    read_number: int,
    output_prefix: FilePath = "",
    gzip: bool = False,
    n_cores: int = 1,
    suffix: str = "",
    **kwargs: Any,
) -> None:
    statement = []
    cat_executable = "zcat"
    split_executable = "split"

    if suffix:
        split_suffix = f"{suffix}_{read_number}.fastq"
    else:
        split_suffix = f"_{read_number}.fastq"

    if ".gz" not in fn:
        cat_executable = "cat"

    if PLATFORM == "darwin":
        if shutil.which("gsplit") is None:
            raise RuntimeError(
                "GNU split is required for unix FASTQ splitting on macOS. "
                "Install coreutils or use --method python."
            )
        split_executable = "gsplit"
        if cat_executable == "zcat":
            cat_executable = "gzcat"

    cmd = (
        f"{cat_executable} {fn} | "
        f"{split_executable} FILTER -l {n_reads * 4} -d "
        f"--additional-suffix={split_suffix} - {output_prefix}_part;"
    )
    if gzip:
        cmd = cmd.replace("FILTER", f"--filter='pigz -p {n_cores} > $FILE.gz'")
    else:
        cmd = cmd.replace("FILTER", "")

    statement.append(cmd)

    logger.info(f"Running: {cmd}")
    subprocess.run(" ".join(statement), shell=True, check=True)


def split_fastq(
    input_files: Sequence[FilePath],
    method: Literal["python", "unix", "seqkit"] = "unix",
    split_type: Literal["n-reads", "n-parts"] = "n-reads",
    output_prefix: FilePath = "split",
    compression_level: int = 5,
    n_reads: int = 1000000,
    n_parts: int = 1,
    suffix: str = "",
    gzip: bool = True,
    n_cores: int = 1,
) -> None:
    """Split FASTQ file(s) into chunks."""

    from capcruncher.api.io import (
        FastqReaderProcess,
        FastqReadFormatterProcess,
        FastqWriterSplitterProcess,
    )

    if split_type == "n-reads" and method == "python":
        readq = SimpleQueue()
        writeq = SimpleQueue()

        paired = len(input_files) > 1

        reader = FastqReaderProcess(
            input_files=input_files,
            outq=readq,
            read_buffer=n_reads,
        )

        formatter = [
            FastqReadFormatterProcess(inq=readq, outq=writeq) for _ in range(1)
        ]

        writer = FastqWriterSplitterProcess(
            inq=writeq,
            output_prefix=output_prefix,
            paired_output=paired,
            n_subprocesses=1,
            gzip=gzip,
            compression_level=compression_level,
        )

        processes = [writer, reader, *formatter]

        for proc in processes:
            proc.start()

        for proc in processes:
            proc.join()
            proc.terminate()

    elif split_type == "n-reads" and method == "unix":
        tasks = []
        n_cores_per_task = (n_cores // 2) if (n_cores // 2) > 1 else 1

        if "," in input_files[0]:
            input_files = [str(fnames).replace(",", " ") for fnames in input_files]

        for ii, fn in enumerate(input_files):
            tasks.append(
                delayed(run_unix_split)(
                    fn,
                    n_reads=n_reads,
                    read_number=ii + 1,
                    gzip=gzip,
                    compression_level=compression_level,
                    output_prefix=output_prefix,
                    n_cores=n_cores_per_task,
                    suffix=suffix,
                )
            )

        Parallel(n_jobs=2 if n_cores > 1 else 1)(tasks)

        for fn in glob.glob(f"{output_prefix}_part*"):
            src = fn
            part_no = int(
                re.match(r"(?:.*)_part(\d+)_.*([1|2])?.fastq(.gz)?", fn).group(1)
            )
            dest = re.sub(r"_part\d+_", f"_part{part_no}_", src)
            os.rename(src, dest)


def digest_fastq(
    fastqs: Sequence[FilePath],
    restriction_site: str,
    mode: str = "pe",
    output_file: Path | str = "out.fastq.gz",
    minimum_slice_length: int = 18,
    statistics: Path | str = "digest.json",
    sample_name: str = "sampleX",
    **kwargs: Any,
) -> Any:
    """Digest FASTQ files and write digestion statistics."""

    from capcruncher_tools.api import digest_fastq as digest_fastq_records

    from capcruncher.utils import get_restriction_site

    logger.info("Digesting FASTQ files")

    if len(fastqs) > 1 and mode == "flashed":
        raise ValueError("Flashed mode can only be used with a single FASTQ file")

    stats = digest_fastq_records(
        fastqs=fastqs,
        restriction_site=get_restriction_site(restriction_site),
        output=output_file,
        read_type=mode.title(),
        sample_name=sample_name,
        minimum_slice_length=minimum_slice_length,
    )

    logger.info("Digestion complete. Generating statistics")
    with open(statistics, "w") as f:
        f.write(stats.model_dump_json())

    return stats


def deduplicate_fastq(
    fastq_1: Sequence[FilePath],
    fastq_2: Sequence[FilePath],
    output_prefix: str | Path = "deduplicated_",
    statistics: str = "deduplication_statistics.json",
    sample_name: str = "sampleX",
    shuffle: bool = False,
    **kwargs: Any,
) -> None:
    """Deduplicate paired FASTQ files and write deduplication statistics."""

    from capcruncher.api.statistics import FastqDeduplicationStatistics
    from capcruncher_tools.api import deduplicate_fastq as deduplicate_fastq_records

    df_stats = deduplicate_fastq_records(
        fastq1=fastq_1,
        fastq2=fastq_2,
        output_prefix=output_prefix,
        sample_name=sample_name,
        shuffle=shuffle,
    )

    dedup_stats = FastqDeduplicationStatistics(
        sample=sample_name,
        total=df_stats.query("stat_type == 'read_pairs_total'")["stat"].values[0],
        duplicates=df_stats.query("stat_type == 'read_pairs_duplicated'")[
            "stat"
        ].values[0],
    )
    with open(statistics, "w") as f:
        f.write(dedup_stats.model_dump_json())

    logger.info("Printing deduplication statistics to stdout")
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(df_vis.to_string(index=False))
