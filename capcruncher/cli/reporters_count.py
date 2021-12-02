import xopen
import os
import tqdm
import logging

from capcruncher.tools.io import (
    CCHDF5ReaderProcess,
    FragmentCountingProcess,
    CCHDF5WriterProcess,
)
from capcruncher.tools.count import get_counts_from_tsv, get_counts_from_tsv_by_batch
from capcruncher.utils import get_categories_from_hdf5_column, get_file_type


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

    input_file_type = get_file_type(reporters) if file_type == "auto" else file_type
    output_file_type = get_file_type(output)

    if output_file_type == "tsv":
        with xopen.xopen(output, mode="wb", threads=4) as writer:

            # Write output file header.
            header = "\t".join(["bin1_id", "bin2_id", "count"]) + "\n"
            writer.write(header.encode())

            if input_file_type == "tsv":

                if low_memory:
                    counts = get_counts_from_tsv_by_batch(
                        reporters=reporters,
                        chunksize=chunksize,
                        remove_capture=remove_capture,
                        remove_exclusions=remove_exclusions,
                        subsample=subsample,
                    )
                else:
                    counts = get_counts_from_tsv(
                        reporters=reporters,
                        remove_capture=remove_capture,
                        remove_exclusions=remove_exclusions,
                        subsample=subsample,
                    )

                for (rf1, rf2), count in counts.items():
                    line = "\t".join([str(rf1), str(rf2), str(count)]) + "\n"
                    writer.write(line.encode())
            else:
                raise NotImplementedError("Currently only tsv -> tsv supported")

    elif output_file_type == "hdf5":

        import multiprocessing

        n_cores = 6
        viewpoints = get_categories_from_hdf5_column(reporters, "slices", "viewpoint")

        viewpoints_queue = multiprocessing.Queue()
        slices_queue = multiprocessing.Queue()
        counts_queue = multiprocessing.Queue()

        reader = CCHDF5ReaderProcess(
            path=reporters, key="slices", inq=viewpoints_queue, outq=slices_queue
        )
        counters = [
            FragmentCountingProcess(
                inq=slices_queue,
                outq=counts_queue,
                remove_capture=remove_capture,
                remove_exclusions=remove_exclusions,
                subsample=subsample,
            )
            for i in range(n_cores)
        ]

        if output_as_cooler:
            writer = CCHDF5WriterProcess(
                inq=counts_queue,
                output_path=output,
                output_format="cooler",
                restriction_fragment_map=fragment_map,
                viewpoint_path=viewpoint_path,
            )
        else:
            writer = CCHDF5WriterProcess(
                inq=counts_queue,
                output_path=output,
                output_format=output_file_type,
                output_key="slices",
                single_file=True,
            )


        processes = [writer, *counters, reader]
        for process in processes:
            process.start()

        for vp in viewpoints:
            viewpoints_queue.put(vp)
        
        viewpoints_queue.put(None)
        reader.join()

        for i in range(n_cores):
            slices_queue.put((None, None))

        for counter in counters:
            counter.join()

        counts_queue.put((None, None))
        writer.join()
