#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: asmith

Simple script to validate run parameters. Verifies that trimming and aligning
occur with user inputs. Test files are provided in PATH_TO_REPO/capture-c/test.
"""

import os
import sys
import yaml
import subprocess as sub
import glob

SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR, SCRIPT_NAME = os.path.split(os.path.dirname(SCRIPT_PATH))

def get_parser(parser=None):
    return parser


def load_yaml(yml):
    with open('capturec_pipeline.yml') as r:
        content = r.read().replace('\t', ' ' * 4)
    return yaml.safe_load(content)

def format_run_cmd(cmd):
    if isinstance(cmd, str):
        return ' '.join([x.strip() for x in cmd.split()])
    elif isinstance(cmd, list):
        return ' '.join([x.strip() for x in cmd])

def test_cmd(cmd, task, **kwargs):
    with open(f'{task}_test.stderr', 'w') as err:
        try:
            sub.check_call(cmd, stderr=err, **kwargs)
            print(f'''\n{task} completed successfuly.\n\nRun parameters: {format_run_cmd(cmd)}''')
            return True
        except Exception as e:
            print(f'Error with {task}:\n\n{e}\n\nSee stderr for further details\n')

def test_trim(parameters, test_dir='test'):
    task = 'Trim_galore'
    cmd = f'''trim_galore
          --cores {parameters["run_options"]["threads"]}
          --paired
          {parameters["trim"]["options"]}
          -o {test_dir}/trim_test/
          {test_dir}/test_1.fastq.gz
          {test_dir}/test_2.fastq.gz
          '''
    return test_cmd(cmd.split(), task)

def test_align(parameters, test_dir='test'):
    task = 'Aligning'
    out_dir = f'{test_dir}/align_test'

    if not os.path.exists(out_dir):
        os.mkdir(out_dir)

    cmd = f'''{parameters["align"]["aligner"]}
              {parameters["align"]["options"]}
              {parameters["align"]["index_flag"]}
              {parameters["align"]["index"]}
              {test_dir}/test_1.fastq.gz
              | samtools view -b > {out_dir}/test.bam
              && samtools sort {out_dir}/test.bam -o {out_dir}/test.sorted.bam -m 2G -@ {parameters["run_options"]["threads"]}
              && mv {out_dir}/test.sorted.bam {out_dir}/test.bam'''

    cmd = ' '.join([x.strip() for x in cmd.split() if not x == 'None'])
    return test_cmd(cmd, task, shell=True)

def main():

    params = load_yaml('capturec_pipeline.yml')
    test_dir = os.path.join(SCRIPT_DIR, 'test')

    tasks = [test_trim, test_align]
    if not all([task(params, test_dir=test_dir) for task in tasks]): # check if all tasks succeed
        pass
    else:
        for err in glob.glob('*.stderr'):
            os.remove(err)

if __name__ == '__main__':
    main()
