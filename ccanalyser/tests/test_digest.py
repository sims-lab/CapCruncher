import pysam
import pytest
import os

from click.testing import CliRunner
from ccanalyser.tools.digest import DigestedRead, DigestedChrom, ReadDigestionProcess
from ccanalyser.utils import PysamFakeEntry
from multiprocessing import SimpleQueue
from multiprocessing import Manager


def test_digest_read_flashed_correct():

    # Correct
    sequence = ''.join(['A'*20, 'G'*20, 'GATC' ,'T'*20, 'C'*20 ])
    read = PysamFakeEntry('read', sequence, ''.join(['A' for _ in sequence]))
    dr_flashed = DigestedRead(read, 'GATC')
    dr_indexes = dr_flashed.get_recognition_site_indexes()
    dr_string = str(dr_flashed)
    dr_string_list = dr_string.split('\n')[:-1]

    assert dr_indexes == [0, 40, 84]
    assert len(dr_string_list) == (2 * 4)

def test_digest_read_too_small():

    # First slice too small
    sequence = ''.join(['A'*5, 'G'*5, 'GATC' ,'T'*20, 'C'*20 ])
    read = PysamFakeEntry('read', sequence, ''.join(['A' for _ in sequence]))
    dr = DigestedRead(read, 'GATC')
    dr_indexes = dr.get_recognition_site_indexes()
    dr_string = str(dr)
    dr_string_list = dr_string.split('\n')[:-1]

    assert len(dr_string_list) == (1 * 4)
    assert len(dr_string_list[1]) == 40

def test_digest_read_no_cutsite():

    #Allowing not digested
    sequence = ''.join(['A'*20, 'G'*20, 'T'*20, 'C'*20 ])
    read = PysamFakeEntry('read', sequence, ''.join(['A' for _ in sequence]))
    dr = DigestedRead(read, 'GATC', allow_undigested=True)
    dr_indexes = dr.get_recognition_site_indexes()
    dr_string = str(dr)
    dr_string_list = dr_string.split('\n')[:-1]

    
    assert dr_indexes == [0, 80]
    assert len(dr_string_list) == (1 * 4)

    #Not allowing undigested
    sequence = ''.join(['A'*20, 'G'*20, 'T'*20, 'C'*20 ])
    read = PysamFakeEntry('read', sequence, ''.join(['A' for _ in sequence]))
    dr = DigestedRead(read, 'GATC', allow_undigested=False)
    dr_indexes = dr.get_recognition_site_indexes()
    dr_string = str(dr)
    dr_string_list = dr_string.split('\n')[:-1]

    
    assert dr_indexes == [0, 80]
    assert len(dr_string_list) == 0

def test_digestion_process_flashed():

    inq = SimpleQueue()
    outq = SimpleQueue()

    #flashed
    dp = ReadDigestionProcess(inq, outq, cutsite='GATC')

    sequence_correct = ''.join(['A'*20, 'G'*20, 'GATC' ,'T'*20, 'C'*20 ])
    sequence_too_short = ''.join(['A'*5, 'G'*5, 'GATC' ,'T'*20, 'C'*20 ])
    sequence_no_cutsite = ''.join(['A'*20, 'G'*20, 'T'*20, 'C'*20 ])
    sequences = [sequence_correct, sequence_too_short, sequence_no_cutsite]

    
    # Correct
    reads_correct = [(PysamFakeEntry(f'read{i}', sequence_correct, ''.join(['A' for _ in sequence_correct])), )
             for i in range(10) ]
    
    # Correct multiple reads
    reads_correct_multiple = [(PysamFakeEntry(f'read{i}', sequence_correct, ''.join(['A' for _ in sequence_correct])),
                     PysamFakeEntry(f'read{i+1}', sequence_correct, ''.join(['A' for _ in sequence_correct])))
             for i in range(10)]
    
    # Short
    reads_short = [(PysamFakeEntry(f'read{i}', sequence_too_short, ''.join(['A' for _ in sequence_too_short])), )
             for i in range(10) ]
    
    # No cutsite
    reads_no_cut = [(PysamFakeEntry(f'read{i}', sequence_no_cutsite, ''.join(['A' for _ in sequence_no_cutsite])), )
             for i in range(10) ]

    
    dp.start()

    # Correct
    inq.put(reads_correct)
    result = outq.get()
    dr_string_list = result.split('\n')[:-1]
    assert len(dr_string_list) == 10 * 4 * 2 # 10 reads, 4 lines per read,  2 slices each

    # Correct multiple
    inq.put(reads_correct_multiple)
    result = outq.get()
    dr_string_list = result.split('\n')[:-1]
    assert len(dr_string_list) == 10 * 4 * 2  * 2 # 10 reads, 4 lines per read,  2 slices each, 2 reads

    # Too short
    inq.put(reads_short)
    result = outq.get()
    dr_string_list = result.split('\n')[:-1]
    assert len(dr_string_list) == 10 * 4 * 1 # 10 reads, 4 lines per read,  1 slices each

    # No cut
    inq.put(reads_no_cut)
    result = outq.get()
    dr_string_list = result.split('\n')[:-1]
    assert len(dr_string_list) == 10 * 4 * 0 # 10 reads, 4 lines per read,  0 slices each

    inq.put('END')
    dp.join()
    dp.terminate()

    
    