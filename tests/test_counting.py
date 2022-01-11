# import os

# import pytest
# from capcruncher.tools.io import (
#     CCHDF5ReaderProcess,
#     CCParquetReaderProcess,
#     CCHDF5WriterProcess,
#     CCCountsWriterProcess,
# )
# import multiprocessing


# @pytest.fixture(scope="module")
# def data_path():
#     fn = os.path.realpath(__file__)
#     dirname = os.path.dirname(fn)
#     data_dir = os.path.join(dirname, "data", "reporter_count")
#     return data_dir

# @pytest.fixture(scope="module")
# def data_parquet(data_path):
#     return os.path.join(data_path, "SAMPLE-A_REP1.parquet")


# @pytest.mark.parametrize(
#     "data,reader,viewpoint,mode,n_records",
#     [   
#         # pytest.param(CCHDF5ReaderProcess, "Slc25A37", 76, id="ok_viewpoint"),
#         # pytest.param(CCHDF5ReaderProcess, "XXXX", 0, id="bad_viewpoint", marks=pytest.mark.xfail),
#         pytest.param(pytest.lazy_fixture("data_parquet"), CCParquetReaderProcess, "Slc25A37", "single", 213, id="ok_viewpoint"),
#         #pytest.param(pytest.lazy_fixture("data_parquet"), CCParquetReaderProcess, "XXXX", 0, id="bad_viewpoint", marks=pytest.mark.xfail),
#     ],
# )
# def test_reader(data, reader, viewpoint,mode, n_records):
    
#     inq = multiprocessing.Queue()
#     outq = multiprocessing.Queue()
#     reader_process = reader(data, inq, outq, selection_mode=mode)

#     reader_process.start()
#     inq.put(viewpoint)
#     (vp, records) = outq.get()
#     inq.put(None)
#     reader_process.join()

#     assert records.shape[0] == n_records




