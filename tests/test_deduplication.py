import glob
import os
import pickle

import pandas as pd
import pytest
import ujson
from capcruncher.cli import fastq_deduplicate
from capcruncher.tools.deduplicate import (
    ReadDeduplicationParserProcess,
    ReadDuplicateRemovalProcess,
)
from capcruncher.tools.io import FastqReaderProcess
from capcruncher.utils import MockFastqRecord, get_file_type, load_dict, save_dict
from pybedtools.bedtool import BedTool
import pysam


@pytest.fixture(scope="module")
def data_path():
    fn = os.path.realpath(__file__)
    dirname = os.path.dirname(fn)
    data_dir = os.path.join(dirname, "data", "fastq_deduplication")
    return data_dir


@pytest.mark.parametrize(
    "fastq_files,outfile,format,n_records_expected",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "parsed.json",
            "json",
            1520,
        ),
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "parsed.pkl",
            "pickle",
            1520,
        ),
    ],
)
def test_fastq_parsing(
    data_path, tmpdir, fastq_files, outfile, format, n_records_expected
):

    infiles = [os.path.join(data_path, fn) for fn in fastq_files]

    outfile_path = os.path.join(tmpdir, outfile)

    fastq_deduplicate.parse(infiles, outfile_path)

    assert os.path.exists(outfile_path)

    with open(outfile_path, "rb") as output:
        if format == "json":
            parse_test = ujson.load(output)
        elif format == "pickle":
            parse_test = pickle.load(output)

    assert len(parse_test) == n_records_expected


@pytest.mark.parametrize(
    "infiles,outfile,n_duplicates_expected",
    [
        (("parsed.json",), "deduplicated.json", 538),
        (("parsed.pickle",), "deduplicated.pickle", 538),
        (("parsed.json", "parsed.pickle"), "deduplicated.pickle", 1520),
    ],
)
def test_fastq_duplicate_identification(
    data_path, tmpdir, infiles, outfile, n_duplicates_expected
):
    infiles_paths = [os.path.join(data_path, fn) for fn in infiles]
    outfile_path = os.path.join(tmpdir, outfile)
    duplicated_ids = fastq_deduplicate.identify(infiles_paths, outfile_path)

    assert os.path.exists(outfile_path)

    outfile_type = get_file_type(outfile)
    output_dict = load_dict(outfile_path, format=outfile_type)

    assert len(duplicated_ids) == n_duplicates_expected


@pytest.mark.parametrize(
    "infiles,duplicates,prefix,n_duplicates_expected",
    [
        (
            ("duplicated_1.fastq.gz", "duplicated_2.fastq.gz"),
            "identified.pkl",
            "out",
            538,
        )
    ],
)
def test_fastq_duplicate_removal(
    data_path, tmpdir, infiles, duplicates, prefix, n_duplicates_expected
):

    infiles_paths = [os.path.join(data_path, fn) for fn in infiles]
    ids = os.path.join(data_path, duplicates)
    out_prefix = os.path.join(tmpdir, prefix)

    stats = fastq_deduplicate.remove(
        infiles_paths,
        duplicated_ids=ids,
        output_prefix=out_prefix,
        stats_prefix=os.path.join(tmpdir, "dedup"),
    )

    stats = stats.set_index("stat_type")

    assert stats.at["reads_removed", "stat"] == n_duplicates_expected

    outfiles = glob.glob(f"{out_prefix}*")
    assert len(outfiles) == len(infiles)

    assert (
        len([r for r in pysam.FastxFile(outfiles[0])]) == stats.at["reads_unique", "stat"])


# pytest.fixture(scope="module")
# def identify_test_data(tmpdir):

#     test_data = [(PysamFakeEntry('1_1', 'AAAA', '\\\\'), PysamFakeEntry('1_2', 'TTTT', '\\\\')),
#                  (PysamFakeEntry('2_1', 'AAAA', '\\\\'), PysamFakeEntry('2_2', 'TTTT', '\\\\')),
#                  (PysamFakeEntry('3_1', 'TTTT', '\\\\'), PysamFakeEntry('3_2', 'AAAA', '\\\\')),
#                 ]


#     save_dict(os.path.join(tmpdir), )


#     test_data_dedup = [(PysamFakeEntry('1_1', 'AAAA', '\\\\'), PysamFakeEntry('1_2', 'TTTT', '\\\\')),
#                        (PysamFakeEntry('3_1', 'TTTT', '\\\\'), PysamFakeEntry('3_2', 'AAAA', '\\\\')),
#                 ]


# test_json_path = os.path.join(dir_test, 'test','dup_parse.json')
# test_duplicates_path = os.path.join(dir_test, 'test','duplicates.json')


# def test_fastq_parsing():

#     inq = SimpleQueue()
#     outq = SimpleQueue()

#     rdpp = ReadDeduplicationParserProcess(inq=inq,
#                                           outq=outq,
#                                           save_hashed_dict_path=test_json_path)

#     rdpp.start()

#     inq.put(test_data)
#     inq.put('END')
#     _ = outq.get()

#     rdpp.join()
#     rdpp.terminate()

#     with open(test_json_path) as r:
#         result = ujson.load(r)


#     hash_function = functools.partial(xxhash.xxh64_intdigest, seed=42)
#     result_expected = {str(hash_function(r1.name + r2.name)): hash_function(r1.sequence + r2.sequence) for r1, r2 in test_data}


#     assert len(result) == len(test_data)
#     assert result == result_expected
#     assert len({x for x in result.values()}) == 2

# def test_fastq_identify():

#     runner = CliRunner()
#     result = runner.invoke(
#             cli,
#             [
#                 "fastq",
#                 "deduplicate",
#                 "identify",
#                 test_json_path,
#                 "-o",
#                 test_duplicates_path,
#             ],
#         )


#     with open(test_duplicates_path) as r:
#         result = ujson.load(r)


#     # Tests that the function identifies the duplicates correctly. Currently should be only one duplicate
#     assert len(result) == 1


# def test_fastq_removal():

#     inq = SimpleQueue()
#     outq = SimpleQueue()

#     with open(test_duplicates_path) as r:
#         duplicates = ujson.load(r)
#         duplicates_set = {int(k) for k in duplicates}


#     rdrp = ReadDuplicateRemovalProcess(inq=inq,
#                                        outq=outq,
#                                        duplicated_ids=duplicates_set,
#                                        )

#     rdrp.start()

#     inq.put(test_data)
#     inq.put('END')
#     result = outq.get()

#     rdrp.join()
#     rdrp.terminate()

#     # Correct number of duplicates removed
#     assert len(result)  == len(test_data_dedup)

#     test_not_duplicated = [(r1.name, r2.name) for r1, r2 in result]
#     expected_not_duplicated = [(r1.name, r2.name) for r1, r2 in test_data_dedup]

#     # Correct duplicate names removed
#     assert test_not_duplicated == expected_not_duplicated
