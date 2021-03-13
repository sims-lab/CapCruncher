import os
from typing import Tuple

import click
from ccanalyser.cli.cli_fastq import cli

from multiprocessing import Queue, SimpleQueue
from ccanalyser.tools.digest import ReadDigestionProcess
from ccanalyser.tools.io import FastqReaderProcess, FastqWriterProcess
from ccanalyser.tools.statistics import DigestionStatCollector, DigestionStatistics
from ccanalyser.utils import get_re_site


@cli.command()
@click.argument("input_fastq", nargs=-1, required=True)
@click.option("-r", "--restriction_enzyme", help="Restriction enzyme name or sequence", required=True)
@click.option(
    "-m",
    "--mode",
    help="Digestion mode. Combined (Flashed) or non-combined (PE) read pairs.",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
    required=True,
)
@click.option("-o", "--output_file", default="out.fastq.gz")
@click.option("-p", "--n_cores", default=1, type=click.INT)
@click.option("--minimum_slice_length", default=18, type=click.INT)
@click.option("--keep_cutsite", default=False)
@click.option(
    "--compression_level",
    help="Level of compression for output files",
    default=5,
    type=click.INT,
)
@click.option(
    "--read_buffer",
    help="Number of reads to process before writing to file",
    default=1e5,
    type=click.INT,
)
@click.option("--stats_prefix", help="Output prefix for stats file", default="stats")
@click.option(
    "--sample_name", help="Name of sample e.g. DOX_treated_1", default="sampleX"
)
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

    \b
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
    inputq = SimpleQueue()  # reads are placed into this queue for processing
    writeq = SimpleQueue()  # digested reads are placed into the queue for writing
    statq = Queue()  # stats queue

    # Variables
    cut_site = get_re_site(restriction_enzyme)

    # Checks the submode to see in which mode to run
    if mode == "flashed":

        # If flashed reads, more confident in presence of rf junction.
        # Will not allow undigested reads in this case as probably junk.
        reader = FastqReaderProcess(
            input_files=input_fastq,
            outq=inputq,
            n_subprocesses=n_cores,
            read_buffer=read_buffer,
        )

        digestion_processes = [
            ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=mode,
                allow_undigested=False,  # Prevents outputting undigested reads
                statq=statq,
            )
            for _ in range(n_cores)
        ]

    elif mode == "pe":

        reader = FastqReaderProcess(
            input_files=input_fastq,
            outq=inputq,
            n_subprocesses=n_cores,
            read_buffer=read_buffer,
        )

        digestion_processes = [
            ReadDigestionProcess(
                inq=inputq,
                outq=writeq,
                cutsite=cut_site,
                min_slice_length=minimum_slice_length,
                read_type=mode,
                allow_undigested=True,
                statq=statq,
            )
            for _ in range(n_cores)
        ]

    # Writer process is common to both
    writer = FastqWriterProcess(
        inq=writeq,
        output=output_file,
        n_subprocesses=n_cores,
        compression_level=compression_level,
    )

    # Start all processes
    processes = [writer, reader, *digestion_processes]

    for proc in processes:
        proc.start()

    reader.join()
    writer.join()

    # Collate stats
    print("")
    print("Collating stats")
    collated_stats = DigestionStatCollector(statq, n_cores).get_collated_stats()

    stats = [
        DigestionStatistics(
            sample=sample_name,
            read_type=mode,
            read_number=read_number,
            slices_unfiltered=stats["unfiltered"],
            slices_filtered=stats["filtered"],
        )
        for read_number, stats in collated_stats.items()
    ]

    if len(stats) > 1:  # Need to collate stats from digestion of 2+ files
        for stat in stats[1:]:
            digestion_stats = stats[0] + stat
    else:
        digestion_stats = stats[0]

    digestion_stats.unfiltered_histogram.to_csv(
        f"{stats_prefix}.digestion.unfiltered.histogram.csv", index=False
    )
    digestion_stats.filtered_histogram.to_csv(
        f"{stats_prefix}.digestion.filtered.histogram.csv", index=False
    )
    digestion_stats.slice_summary.to_csv(
        f"{stats_prefix}.digestion.slice.summary.csv", index=False
    )
    digestion_stats.read_summary.to_csv(
        f"{stats_prefix}.digestion.read.summary.csv", index=False
    )

    print(digestion_stats.read_summary)
