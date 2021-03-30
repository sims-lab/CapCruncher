from typing import Dict, TypedDict, Union
import pybedtools
from pybedtools.bedtool import BedTool
import pysam
import pytest
import os
import hashlib
import xxhash
import ujson
import pickle
import operator
from click.testing import CliRunner
import functools
import pandas as pd

from ccanalyser.cli import cli
from ccanalyser.tools.deduplicate import ReadDeduplicationParserProcess, ReadDuplicateRemovalProcess
from multiprocessing import SimpleQueue, Queue
from ccanalyser.utils import PysamFakeEntry

dir_test = os.path.realpath(os.path.dirname(__file__))
dir_package = os.path.dirname(dir_test)
dir_data = os.path.join(dir_package, "data")

if not os.path.exists(os.path.join(dir_test, 'test')):
    os.mkdir(os.path.join(dir_test, 'test'))

if not os.path.exists(os.path.join(dir_test, 'stats')):
    os.mkdir(os.path.join(dir_test, 'stats'))


test_data = [(PysamFakeEntry('1_1', 'AAAA', '\\\\'), PysamFakeEntry('1_2', 'TTTT', '\\\\')),
             (PysamFakeEntry('2_1', 'AAAA', '\\\\'), PysamFakeEntry('2_2', 'TTTT', '\\\\')),
             (PysamFakeEntry('3_1', 'TTTT', '\\\\'), PysamFakeEntry('3_2', 'AAAA', '\\\\')),
            ]

test_data_dedup = [(PysamFakeEntry('2_1', 'AAAA', '\\\\'), PysamFakeEntry('2_2', 'TTTT', '\\\\')),
                   (PysamFakeEntry('3_1', 'TTTT', '\\\\'), PysamFakeEntry('3_2', 'AAAA', '\\\\')),
            ]

test_json_path = os.path.join(dir_test, 'test','dup_parse.json')
test_duplicates_path = os.path.join(dir_test, 'test','duplicates.json')


def test_fastq_parsing():

    inq = SimpleQueue()
    outq = SimpleQueue()

    rdpp = ReadDeduplicationParserProcess(inq=inq, 
                                          outq=outq, 
                                          save_hashed_dict_path=test_json_path)

    rdpp.start()
    
    inq.put(test_data)
    inq.put('END')
    _ = outq.get()

    rdpp.join()
    rdpp.terminate()

    with open(test_json_path) as r:
        result = ujson.load(r)


    hash_function = functools.partial(xxhash.xxh64_intdigest, seed=42)
    result_expected = {str(hash_function(r1.name + r2.name)): hash_function(r1.sequence + r2.sequence) for r1, r2 in test_data}


    assert len(result) == len(test_data)
    assert result == result_expected
    assert len({x for x in result.values()}) == 2 

def test_fastq_identify():
    
    runner = CliRunner()
    result = runner.invoke(
            cli,
            [
                "fastq",
                "deduplicate",
                "identify",
                test_json_path,
                "-o",
                test_duplicates_path,
            ],
        )

    
    with open(test_duplicates_path) as r:
        result = ujson.load(r)
    

    # Tests that the function identifies the duplicates correctly. Currently should be only one duplicate
    assert len(result) == 1


    
def test_fastq_removal():

    inq = SimpleQueue()
    outq = SimpleQueue()

    with open(test_duplicates_path) as r:
        duplicates = ujson.load(r)
        duplicates_set = set(duplicates)


    rdrp = ReadDuplicateRemovalProcess(inq=inq, 
                                       outq=outq, 
                                       duplicated_ids=duplicates_set,
                                       )

    rdrp.start()
    
    inq.put(test_data)
    inq.put('END')
    result = outq.get()

    rdrp.join()
    rdrp.terminate()


    # Correct number of duplicates removed
    assert len(result)  == len(test_data_dedup)
    
    test_not_duplicated = [(r1.name, r2.name) for r1, r2 in result]
    expected_not_duplicated = [(r1.name, r2.name) for r1, r2 in test_data_dedup]

    # Correct duplicate names removed
    assert test_not_duplicated == expected_not_duplicated



    


    




