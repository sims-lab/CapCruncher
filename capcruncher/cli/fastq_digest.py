import os
from typing import Tuple
from multiprocessing import Queue
from capcruncher.tools.digest import ReadDigestionProcess
from capcruncher.tools.io import FastqReaderProcess, FastqWriterProcess
from capcruncher.utils import get_re_site
import pandas as pd


def collate_statistics(statq: Queue, 
                       n_subprocesses: int) -> pd.DataFrame:
    """Collates digestion statistics from supplied statistics queue.

    Args:
        statq (Queue): Queue to use for collating statistics. Final item(s) must be 'END'
        n_subprocesses (int): Number of digestion processes used. Required to know when to stop aquiring
                              data from the queue.

    Returns:
        pd.DataFrame: Digestion statistics in histogram format. 
                      Columns: 'read_type', 'read_number', 'unfiltered/filtered', 'n_slices', 'n_reads'
    """

    n_subprocesses_terminated = 0
    dframes = []

    while n_subprocesses_terminated < n_subprocesses:

        stats = statq.get()
           
        if stats == 'END':
            n_subprocesses_terminated += 1
            continue
    
        else:
            df = pd.DataFrame(stats)
            df_melt = df.melt(id_vars=['read_type', 'read_number'], var_name='unfiltered/filtered', value_name='n_slices')
            df_hist = (df_melt.groupby(['read_type', 'read_number', 'unfiltered/filtered', 'n_slices'])
                              .size()
                              .to_frame('n_reads')
                              .reset_index()
                              .sort_values('n_slices'))

            dframes.append(df_hist)

    
    return (pd.concat(dframes)
              .groupby(['read_type', 'read_number', 'unfiltered/filtered', 'n_slices'])
              .sum()
              .reset_index()
              .sort_values('n_slices')
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
    statq = Queue()  # stats queue

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


    # Collate statistics
    df_stats = collate_statistics(statq=statq, n_subprocesses=n_cores)
    df_stats['sample'] = sample_name 

    # Ensure that all processes have rejoined the main thread.
    reader.join()
    writer.join()
    
    # Ensure that the queues are joined to the main process
    inputq.close()
    writeq.close()
    statq.close()

    # Kill the digestion processes
    for proc in digestion_processes:
        proc.terminate()


    # Unfiltered histogram
    (df_stats.query('`unfiltered/filtered` == "unfiltered"')
            .to_csv(f"{stats_prefix}.digestion.unfiltered.histogram.csv", index=False)
    )

    # Filtered histogram
    (df_stats.query('`unfiltered/filtered` == "filtered"')
            .to_csv(f"{stats_prefix}.digestion.filtered.histogram.csv", index=False)
    )

    # Read summary - reads that pass the filter
    (df_stats.query('n_slices > 0')
             .groupby(['sample', 'read_type', 'read_number', 'unfiltered/filtered'])
             ['n_reads']
             .sum()
             .reset_index()
             .assign(stage='digestion')
             .rename(columns={'unfiltered/filtered': 'stat_type', 'n_reads': 'stat'})
             .to_csv(f"{stats_prefix}.digestion.read.summary.csv", index=False)
    )

