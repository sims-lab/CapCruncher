from builtins import breakpoint
import multiprocessing
import os
from typing import Tuple
from multiprocessing import Queue
from capcruncher.tools.digest import ReadDigestionProcess
from capcruncher.tools.io import FastqReaderProcess, FastqWriterProcess
from capcruncher.utils import get_re_site
import pandas as pd

def digest(
    input_fastq: Tuple,
    restriction_enzyme: str,
    mode: str = "pe",
    output_file: os.PathLike = "out.fastq.gz",
    minimum_slice_length: int = 18,
    compression_level: int = 5,
    n_cores: int = 1,
    read_buffer: int = 100000,
    stats_prefix: os.PathLike = "",
    keep_cutsite: bool = False,
    sample_name: str = "",
):
    """
    Performs in silico digestion of one or a pair of fastq files.

    \f
    Args:
     input_fastq (Tuple): Input fastq files to process
     restriction_enzyme (str): Restriction enzyme name or site to use for digestion.
     mode (str, optional): Digest combined(flashed) or non-combined(pe).
                           Undigested pe reads are output but flashed are not written. Defaults to "pe".
     output_file (os.PathLike, optional): Output fastq file path. Defaults to "out.fastq.gz".
     minimum_slice_length (int, optional): Minimum allowed length for in silico digested reads. Defaults to 18.
     compression_level (int, optional): Compression level for gzip output (1-9). Defaults to 5.
     n_cores (int, optional): Number of digestion processes to use. Defaults to 1.
     read_buffer (int, optional): Number of reads to process before writing to file. Defaults to 100000.
     stats_prefix (os.PathLike, optional): Output prefix for stats file. Defaults to "".
     keep_cutsite (bool, optional): Determines if cutsite is removed from the output. Defaults to False.
     sample_name (str, optional): Name of sample processed eg. DOX-treated_1. Defaults to ''.
    """

    # Set up multiprocessing variables
    inputq = Queue()  # reads are placed into this queue for processing
    writeq = Queue()  # digested reads are placed into the queue for writing

    cut_site = get_re_site(restriction_enzyme)
    n_digestion_processes = (n_cores - 2) if n_cores > 2 else 1

    reader = FastqReaderProcess(
        input_files=input_fastq,
        outq=inputq,
        read_buffer=read_buffer,
    )

    digestion_processes = list()
    for _ in range(n_digestion_processes):
        stats_pipe = multiprocessing.Pipe()

        if mode == "flashed" and len(input_fastq) == 1:
            digestion_process = ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=mode,
                allow_undigested=False,  # Prevents outputting undigested reads
                stats_pipe=stats_pipe[0],
            )
        elif mode == "pe" and len(input_fastq) == 2:
            digestion_process = ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=mode,
                allow_undigested=True,
                stats_pipe=stats_pipe[0],
            )
        
        else:
            raise ValueError("Wrong number of files specified. Flashed == 1, PE == 2")

        digestion_processes.append((stats_pipe, digestion_process))

    # Writer process is common to both
    writer = FastqWriterProcess(
        inq=writeq,
        output=output_file,
        compression_level=compression_level,
    )

    # Start all processes
    reader.start()
    writer.start()

    for (_, process) in digestion_processes:
        process.start()

    reader.join()

    for _ in range(n_cores - 2):
        inputq.put_nowait(None)

    stats = list()
    for (stats_pipe, digestion_process) in digestion_processes:
        stats.append(stats_pipe[1].recv())
        digestion_process.join()

    writeq.put_nowait(None)

    writer.join()

    # Collate statistics
    df_hist = (
        pd.concat([df for df in stats if df is not None])
        .groupby(["read_type", "read_number", "unfiltered", "filtered"])
        .sum()
        .reset_index()
    )

    # Melt dataframe for manipulation
    df_hist = (
        df_hist.reset_index()
        .melt(
            id_vars=["index", "read_type", "read_number"],
            value_vars=["filtered", "unfiltered"],
            var_name="filtered",
            value_name="n_slices",
        )
        .replace("filtered", True)
        .replace("unfiltered", False)
        .set_index("index")
        .join(df_hist["count"])
        .groupby(["read_type", "read_number", "filtered", "n_slices"])
        .sum()
        .reset_index()
        .assign(sample=sample_name)
    )

    # Unfiltered histogram
    df_hist_unfilt = (df_hist.query("filtered == False")).sort_values(
        ["n_slices", "count"]
    )

    # Filtered histogram
    df_hist_filt = (df_hist.query("filtered == True and n_slices > 0")).sort_values(
        ["n_slices", "count"]
    )

    #breakpoint()

    # Read summary - reads that pass the filter
    df_stats = (
        df_hist
        .query("(n_slices >= 1) or (filtered == False) or (read_type == 'pe')")
        .query("read_number < 2")
        .groupby(["sample", "read_type", "read_number", "filtered"])["count"]
        .sum()
        .reset_index()
        .assign(
            stage="digestion",
            filtered=lambda df: df["filtered"]
            .replace(True, "filtered")
            .replace(False, "unfiltered"),
        )
        .rename(columns={"filtered": "stat_type", "count": "stat"})
    )

    df_hist_unfilt.to_csv(
        f"{stats_prefix}.digestion.unfiltered.histogram.csv", index=False
    )
    df_hist_filt.to_csv(f"{stats_prefix}.digestion.filtered.histogram.csv", index=False)
    df_stats.to_csv(f"{stats_prefix}.digestion.read.summary.csv", index=False)

    return df_stats
