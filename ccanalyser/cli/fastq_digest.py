from multiprocessing import Queue, SimpleQueue

import click
import numpy as np
import pandas as pd
from pysam import FastxFile
from xopen import xopen


from ccanalyser.cli import cli
from ccanalyser.tools.digest import ReadDigestionProcess
from ccanalyser.tools.io import FastqReaderProcess, FastqWriterProcess
from ccanalyser.tools.statistics import (DigestionStatCollector,
                                         DigestionStatistics)
from ccanalyser.utils import get_re_site



@cli.command()
@click.argument("input_fastq", nargs=-1)
@click.option("-r", "--restriction_enzyme", help="Restriction enzyme name or sequence")
@click.option(
    "-m",
    "--mode",
    help="Digestion mode",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
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
@click.option("--stats_prefix", help="Output prefix for stats file")
@click.option("--sample_name", help="Name of sample e.g. DOX_treated_1")
def fastq_digest(
    input_fastq,
    restriction_enzyme=None,
    mode="pe",
    output_file="out.fastq.gz",
    minimum_slice_length=18,
    compression_level=5,
    n_cores=1,
    read_buffer=100000,
    stats_prefix="",
    keep_cutsite=False,
    sample_name='',
):

    '''Performs in silico digestion of an interleaved fastq file'''


    #TODO: Enable keeping the cutsite

    # Set up multiprocessing variables
    inputq = SimpleQueue()  # reads are placed into this queue for processing
    writeq = SimpleQueue()  # digested reads are placed into the queue for writing
    statq = Queue()  # stats queue

    # Variables
    cut_site = get_re_site(restriction_enzyme)

    if mode == "flashed":  # Checks the submode to see in which mode to run

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
                allow_undigested=False,
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
    #sample_name = re.match(r".*/(.*)(_part\d+.).*", input_fastq[0]).group(1)

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
