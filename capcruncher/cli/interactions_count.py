from concurrent.futures import thread
import xopen
import os
import logging
import tempfile
import glob
import more_itertools
import multiprocessing

import pandas as pd
from capcruncher.api.io import (
    CCParquetReaderProcess,
    FragmentCountingProcess,
    CCCountsWriterProcess,
)
from capcruncher.api.storage import merge_coolers

def get_number_of_reader_threads(n_cores):
    threads = (n_cores // 2)

    if 1 < threads < 4:
        threads = threads
    elif threads < 1:
        threads = 1
    else:
        threads = 4
    
    return threads




def count(
    reporters: os.PathLike,
    output: os.PathLike = "counts.tsv",
    file_type: str = "auto",
    remove_exclusions: bool = False,
    remove_capture: bool = False,
    subsample: int = 0,
    low_memory: bool = False,
    chunksize: int = 2e6,
    output_as_cooler: bool = False,
    fragment_map: os.PathLike = None,
    viewpoint_path: os.PathLike = None,
    n_cores: int = 1,
):
    """
    Determines the number of captured restriction fragment interactions genome wide.

    Parses a reporter slices tsv and counts the number of unique restriction fragment
    interaction combinations that occur within each fragment.

    Options to ignore unwanted counts e.g. excluded regions or capture fragments are provided.
    In addition the number of reporter fragments can be subsampled if required.

    \f
    Args:
     reporters (os.PathLike): Reporter tsv file path.
     output (os.PathLike, optional): Output file path for interaction counts tsv. Defaults to 'counts.tsv'.
     remove_exclusions (bool, optional): Removes regions marked as capture proximity exclusions before counting. Defaults to False.
     remove_capture (bool, optional): Removes all capture fragments before counting. Defaults to False.
     subsample (int, optional): Subsamples the fragments by the specified fraction. Defaults to 0 i.e. No subsampling.
    """

    import pyarrow
    import dask.dataframe as dd

    logging.info(f"Examining viewpoints from parquet file: {reporters}")
    # Unsure of the best way to do this. Will just load the first partion vp column and extract
    ddf = dd.read_parquet(
        reporters,
        columns=[
            "viewpoint",
        ],
        engine="pyarrow"
    )

    viewpoints_col = ddf["viewpoint"].cat.as_known()
    viewpoint_sizes = viewpoints_col.value_counts().compute()
    viewpoints = list(viewpoints_col.cat.categories)

    logging.info(f"Number of viewpoints: {len(viewpoints)}")
    logging.info(f"Number of slices per viewpoint: {viewpoint_sizes.to_dict()}")

    MAX_SLICE_NUMBER = 1e6

    if (
        viewpoint_sizes[viewpoint_sizes > MAX_SLICE_NUMBER].shape[0] == 0
    ):  # Less than MAX_SLICE_NUMBER slices, read in batch
        mode = "batch"
    else:
        # More than MAX_SLICE_NUMBER slices, read by partition
        mode = "partition"

    # Multiprocessing queues
    viewpoints_queue = multiprocessing.Queue()
    slices_queue = multiprocessing.Queue()
    counts_queue = multiprocessing.Queue()
    

    reader = CCParquetReaderProcess(
        path=reporters,
        inq=viewpoints_queue,
        outq=slices_queue,
        selection_mode=mode,
    )

    # Multiprocessing set-up
    n_reading_threads = get_number_of_reader_threads(n_cores=n_cores)
    pyarrow.set_cpu_count(n_reading_threads)

    n_worker_processes = ((n_cores - n_reading_threads) // 2)
    n_counting_processes = n_worker_processes if n_worker_processes > 1 else 1
    n_writing_processes =  n_counting_processes

    logging.info(f"Starting {n_reading_threads} reader threads")

    counters = [
        FragmentCountingProcess(
            inq=slices_queue,
            outq=counts_queue,
            remove_capture=remove_capture,
            remove_exclusions=remove_exclusions,
            subsample=subsample,
        )
        for _ in range(n_counting_processes)
    ]

    tmpdir = tempfile.TemporaryDirectory()
    writers = [
        CCCountsWriterProcess(
            inq=counts_queue,
            output_format="cooler" if output_as_cooler else file_type,
            restriction_fragment_map=fragment_map,
            viewpoint_path=viewpoint_path,
            tmpdir=tmpdir.name,
        )
        for i in range(n_writing_processes)
    ]

    # Start all processes
    processes = [*writers, *counters, reader]
    for process in processes:
        process.start()


    for vp_chunk in more_itertools.chunked(viewpoints, 10):
        viewpoints_queue.put(vp_chunk)

    viewpoints_queue.put(None) # Sentinel value to signal end of viewpoints
    reader.join()

    # End the counting inqueue
    for _ in range(n_cores):
        slices_queue.put((None, None))

    # Join the counters
    for counter in counters:
        counter.join()

    # End the counting queue
    for _ in range(n_cores):
        counts_queue.put((None, None))

    # Join the writers
    for writer in writers:
        writer.join()

    # Merge the output files together
    # TODO: Allow other than cooler outputs
    logging.info(f"Making final cooler at {output}")
    output_files = glob.glob(os.path.join(tmpdir.name, "*.hdf5"))
    merge_coolers(output_files, output=output)

    tmpdir.cleanup()
