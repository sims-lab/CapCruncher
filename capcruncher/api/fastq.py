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
from pydantic import BaseModel, Field, PositiveInt, field_validator, model_validator

PLATFORM = sys.platform


def _as_existing_paths(paths: Sequence[Path | str]) -> tuple[Path, ...]:
    normalised_paths = tuple(Path(path) for path in paths)
    missing_paths = [path for path in normalised_paths if not path.exists()]
    if missing_paths:
        missing = ", ".join(str(path) for path in missing_paths)
        raise ValueError(f"Input path(s) do not exist: {missing}")
    return normalised_paths


class FastqSplitOptions(BaseModel):
    """Validated options for FASTQ splitting."""

    input_files: tuple[Path, ...]
    method: Literal["python", "unix", "seqkit"] = "unix"
    split_type: Literal["n-reads", "n-parts"] = "n-reads"
    output_prefix: Path = Path("split")
    compression_level: int = Field(default=5, ge=0, le=9)
    n_reads: PositiveInt = 1_000_000
    n_parts: PositiveInt = 1
    suffix: str = ""
    gzip: bool = True
    n_cores: PositiveInt = 1

    @field_validator("input_files", mode="before")
    @classmethod
    def validate_input_files(cls, value: Sequence[Path | str]) -> tuple[Path, ...]:
        return _as_existing_paths(value)

    @field_validator("input_files")
    @classmethod
    def validate_fastq_count(cls, value: tuple[Path, ...]) -> tuple[Path, ...]:
        if not value:
            raise ValueError("At least one FASTQ file is required.")
        if len(value) > 2:
            raise ValueError("FASTQ splitting accepts one file or one read pair.")
        return value


class FastqDigestOptions(BaseModel):
    """Validated options for FASTQ digestion."""

    fastqs: tuple[Path, ...]
    restriction_site: str = Field(min_length=1)
    mode: Literal["flashed", "pe"] = "pe"
    output_file: Path = Path("out.fastq.gz")
    minimum_slice_length: PositiveInt = 18
    statistics: Path = Path("digest.json")
    sample_name: str = Field(default="sampleX", min_length=1)

    @field_validator("fastqs", mode="before")
    @classmethod
    def validate_fastqs(cls, value: Sequence[Path | str]) -> tuple[Path, ...]:
        return _as_existing_paths(value)

    @model_validator(mode="after")
    def validate_mode_file_count(self) -> "FastqDigestOptions":
        if self.mode == "flashed" and len(self.fastqs) != 1:
            raise ValueError("Flashed mode requires exactly one FASTQ file.")
        if self.mode == "pe" and len(self.fastqs) != 2:
            raise ValueError("PE mode requires exactly two FASTQ files.")
        return self


class FastqDeduplicationOptions(BaseModel):
    """Validated options for paired FASTQ deduplication."""

    fastq_1: tuple[Path, ...]
    fastq_2: tuple[Path, ...]
    output_prefix: Path = Path("deduplicated_")
    statistics: Path = Path("deduplication_statistics.json")
    sample_name: str = Field(default="sampleX", min_length=1)
    shuffle: bool = False

    @field_validator("fastq_1", "fastq_2", mode="before")
    @classmethod
    def validate_fastqs(cls, value: Sequence[Path | str]) -> tuple[Path, ...]:
        return _as_existing_paths(value)

    @model_validator(mode="after")
    def validate_read_pairs(self) -> "FastqDeduplicationOptions":
        if not self.fastq_1 or not self.fastq_2:
            raise ValueError("Both FASTQ read lists are required.")
        if len(self.fastq_1) != len(self.fastq_2):
            raise ValueError("FASTQ read lists must contain the same number of files.")
        return self


def run_unix_split(
    fn: Path,
    n_reads: int,
    read_number: int,
    output_prefix: Path = Path(),
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

    if ".gz" not in str(fn):
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
    input_files: Sequence[Path | str],
    method: Literal["python", "unix", "seqkit"] = "unix",
    split_type: Literal["n-reads", "n-parts"] = "n-reads",
    output_prefix: Path | str = Path("split"),
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

    options = FastqSplitOptions(
        input_files=input_files,
        method=method,
        split_type=split_type,
        output_prefix=Path(output_prefix),
        compression_level=compression_level,
        n_reads=n_reads,
        n_parts=n_parts,
        suffix=suffix,
        gzip=gzip,
        n_cores=n_cores,
    )
    input_files = options.input_files
    method = options.method
    split_type = options.split_type
    output_prefix = options.output_prefix
    compression_level = options.compression_level
    n_reads = options.n_reads
    gzip = options.gzip
    n_cores = options.n_cores
    suffix = options.suffix

    if split_type == "n-reads" and method == "python":
        readq = SimpleQueue()
        writeq = SimpleQueue()

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
            paired_output=len(input_files) > 1,
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

        if "," in str(input_files[0]):
            input_files = tuple(
                Path(str(fnames).replace(",", " ")) for fnames in input_files
            )

        for read_number, fn in enumerate(input_files, start=1):
            tasks.append(
                delayed(run_unix_split)(
                    fn,
                    n_reads=n_reads,
                    read_number=read_number,
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
    fastqs: Sequence[Path | str],
    restriction_site: str,
    mode: str = "pe",
    output_file: Path | str = Path("out.fastq.gz"),
    minimum_slice_length: int = 18,
    statistics: Path | str = Path("digest.json"),
    sample_name: str = "sampleX",
    **kwargs: Any,
) -> Any:
    """Digest FASTQ files and write digestion statistics."""

    from capcruncher_tools.api import digest_fastq as digest_fastq_records

    from capcruncher.utils import get_restriction_site

    options = FastqDigestOptions(
        fastqs=fastqs,
        restriction_site=restriction_site,
        mode=mode,
        output_file=Path(output_file),
        minimum_slice_length=minimum_slice_length,
        statistics=Path(statistics),
        sample_name=sample_name,
    )

    logger.info("Digesting FASTQ files")

    stats = digest_fastq_records(
        fastqs=tuple(str(fastq) for fastq in options.fastqs),
        restriction_site=get_restriction_site(options.restriction_site),
        output=str(options.output_file),
        read_type=options.mode.title(),
        sample_name=options.sample_name,
        minimum_slice_length=options.minimum_slice_length,
    )

    logger.info("Digestion complete. Generating statistics")
    with open(options.statistics, "w") as f:
        f.write(stats.model_dump_json())

    return stats


def deduplicate_fastq(
    fastq_1: Sequence[Path | str],
    fastq_2: Sequence[Path | str],
    output_prefix: Path | str = Path("deduplicated_"),
    statistics: Path | str = Path("deduplication_statistics.json"),
    sample_name: str = "sampleX",
    shuffle: bool = False,
    **kwargs: Any,
) -> None:
    """Deduplicate paired FASTQ files and write deduplication statistics."""

    from capcruncher.api.statistics import FastqDeduplicationStatistics
    from capcruncher_tools.api import deduplicate_fastq as deduplicate_fastq_records

    output_prefix_for_tools = os.fspath(output_prefix)
    options = FastqDeduplicationOptions(
        fastq_1=fastq_1,
        fastq_2=fastq_2,
        output_prefix=Path(output_prefix),
        statistics=Path(statistics),
        sample_name=sample_name,
        shuffle=shuffle,
    )

    df_stats = deduplicate_fastq_records(
        fastq1=tuple(str(fastq) for fastq in options.fastq_1),
        fastq2=tuple(str(fastq) for fastq in options.fastq_2),
        output_prefix=output_prefix_for_tools,
        sample_name=options.sample_name,
        shuffle=options.shuffle,
    )

    dedup_stats = FastqDeduplicationStatistics(
        sample=options.sample_name,
        total=df_stats.query("stat_type == 'read_pairs_total'")["stat"].values[0],
        duplicates=df_stats.query("stat_type == 'read_pairs_duplicated'")[
            "stat"
        ].values[0],
    )
    with open(options.statistics, "w") as f:
        f.write(dedup_stats.model_dump_json())

    logger.info("Printing deduplication statistics to stdout")
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(df_vis.to_string(index=False))
