#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  8 15:45:09 2020

@author: asmith

Script splits a fastq into specified chunks
"""

import itertools
import logging
from multiprocessing import SimpleQueue
from typing import Tuple
import subprocess
import glob
import os
import re
from joblib import Parallel, delayed


def run_unix_split(fn: os.PathLike, 
                   n_reads: int, 
                   read_number: int,
                   output_prefix: os.PathLike = '',
                   gzip: bool = False,
                   compression_level: int = 5,
                   n_cores=1):

    statement = []    

    cmd = f"""zcat {fn} | split FILTER -l {n_reads * 4} -d --additional-suffix=_{read_number}.fastq - {output_prefix}_part;"""
    
    if not ".gz" in fn:
        cmd = cmd.replace("zcat", "cat")
    
    if gzip:
        cmd = cmd.replace("FILTER", f"--filter='pigz -p {n_cores} > $FILE.gz'")
    else:
        cmd = cmd.replace("FILTER", "")
    
    statement.append(cmd)

    logging.info(f"Running: {cmd}")
    subprocess.run(' '.join(statement), shell=True)


def split(
    input_files: Tuple,
    method: str = "unix",
    output_prefix: os.PathLike = "split",
    compression_level: int = 5,
    n_reads: int = 1000000,
    gzip: bool = True,
    n_cores: int = 1,
):
    """ 
    Splits fastq file(s) into equal chunks of n reads.

    Will now need "," between files of the same read.

    \f
    Args:
     input_files (Tuple): Input fastq files to process.
     method (str, optional): Python or unix method (faster but not guarenteed to mantain read pairings) to split the fastq files. Defaults to "unix".
     output_prefix (os.PathLike, optional): Output prefix for split fastq files. Defaults to "split".
     compression_level (int, optional): Compression level for gzipped output. Defaults to 5.
     n_reads (int, optional): Number of reads to split the input fastq files into. Defaults to 1000000.
     gzip (bool, optional): Gzip compress output files if True. Defaults to True.
    
    """

    from capcruncher.tools.io import (
    FastqReaderProcess,
    FastqWriterSplitterProcess,
    FastqReadFormatterProcess,)


    if method == "python":
        readq = SimpleQueue()
        writeq = SimpleQueue()

        paired = True if len(input_files) > 1 else False

        reader = FastqReaderProcess(
            input_files=input_files,
            outq=readq,
            read_buffer=n_reads,
            n_subprocesses=1,
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

    elif method == "unix": # Using unix split to perform the splitting

        tasks = []
        n_cores_per_task = (n_cores // 2) if (n_cores // 2) > 1 else 1

        if "," in input_files[0]: # Allows for specifying multiple files
            input_files = [fnames.replace(",", " ") for fnames in input_files]

        for ii, fn in enumerate(input_files):
            t = delayed(run_unix_split)(fn, 
                                        n_reads=n_reads,
                                        read_number=ii+1,
                                        gzip=gzip,
                                        compression_level=compression_level,
                                        output_prefix=output_prefix,
                                        n_cores=n_cores_per_task)
            
            tasks.append(t)
        
        # Run splitting
        Parallel(n_jobs=2 if n_cores > 1 else 1)(tasks)


        # The suffixes are in the format 00, 01, 02 etc need to replace with int
        for fn in glob.glob(f"{output_prefix}_part*"):
            src = fn
            part_no = int(re.match(r"(?:.*)_part(\d+)_[1|2].*", fn).group(1))
            dest = re.sub(r"_part\d+_", f"_part{part_no}_", src)
            os.rename(src, dest)
