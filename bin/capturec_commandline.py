import argparse
import os
import sys
from functools import wraps

# Make sure the script can find capturec scripts
SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(PACKAGE_DIR)

def add_subcommand(func,
                   parser,
                   subcommand,
                   help=None):

    subparser = parser.add_parser(subcommand, help=help)

    return func(subparser)

def get_args_deduplicate_fastq(parser=None):

    if not parser:
        parser = argparse.ArgumentParser()

    parser.add_argument('-1', '--fq1', help='fastq file to parse containing read 1')
    parser.add_argument('-2', '--fq2', help='fastq file to parse containing read 2')
    parser.add_argument('--out1', help='fastq file to parse containing read 1',
                       default='out1.fq.gz')
    parser.add_argument('--out2', help='fastq file to parse containing read 2',
                       default='out2.fq.gz')
    parser.add_argument('--stats_file', help='name of deduplication statistics file',
                   default=sys.stdout)
    parser.add_argument('-c', '--compression_level',
                        help='Level of compression (1-9 with 9 being the highest)',
                        type=int, default=5)
    parser.add_argument('--read_buffer',
                   help='defines the number of reads processed before writing to file',
                   default=10000, type=int)

    return parser


def add_digest_shared_options(parser):
    parser.add_argument('-o', '--output_file', help='output file name',
                        default='digested.fastq.gz')

    enzyme_group = parser.add_mutually_exclusive_group(required=True)
    enzyme_group.add_argument('-r', '--restriction_enzyme', help='Name of restriction enzyme',
                              default='DpnII')
    enzyme_group.add_argument('-s', '--cut_sequence',
                               help='Sequence of restriction site')

    parser.add_argument('-m', '--minimum_slice_length', help='Shortest length for a slice to be output',
                        default=20, type=int)
    parser.add_argument('--stats_file',
                        help='stats_file_prefix', default='stats.log')
    parser.add_argument('-c', '--compression_level', help='Level of gzip compression (1-9 with 9 being the most compressed/slowest)',
                        default=6, type=int)
    parser.add_argument('--keep_cutsite', help='Determines if cutsite is stripped from the start of each slice',
                        action='store_true', default=False)
    parser.add_argument(
        '--buffer', help='Number of reads to process before writing output', default=10000, type=int)

    parser.add_argument(
        '-p', '--n_digestion_processes',help='Number of digestion processes to spawn', default=1, type=int)

    return parser

def get_args_digest_fastq_flashed(parser=None):

    if not parser:
        parser = argparse.ArgumentParser(prog='digest_fastq_flashed')

    parser.add_argument('-i', '--input_fastq',
                        help='fastq file to parse',
                        required=True)
    parser = add_digest_shared_options(parser)

    return parser

def get_args_digest_fastq_unflashed(parser=None):

    if not parser:
        parser = argparse.ArgumentParser(prog='digest_fastq_unflashed')

    parser.add_argument('-1', '--fq1',
                        help='fastq file containing read 1 to parse',
                        required=True)
    parser.add_argument('-2', '--fq2',
                                  help='fastq file containing read 2 to parse',
                                  required=True)

    parser = add_digest_shared_options(parser)

    return parser


def collect_arguments():
    parser = argparse.ArgumentParser(prog='capturec',
                                     description='Collection of scripts to process NG-Capture-C data.'
                                    )
    subparser = parser.add_subparsers(dest='subcommand')

    sub_commands = {'dedup_fastq': {'func': get_args_deduplicate_fastq,
                                    'help': 'Removes PCR duplicates from fastq files'},

                    'digest_fastq_flashed': {'func': get_args_digest_fastq_flashed,
                                             'help': 'Digestion of fastq file in silico'},

                    'digest_fastq_unflashed': {'func': get_args_digest_fastq_unflashed,
                                               'help': 'Digestion of fastq file in silico'}



                   }

    for command, details in sub_commands.items():
        subparser = add_subcommand(func=details['func'],
                                   parser=subparser,
                                   subcommand=command,
                                   help=details['help'])


    return parser.parse_args()

def run():

    args = collect_arguments()

    if args.subcommand == 'dedup_fastq':
        from capturec import deduplicate_fastq
        del args.subcommand
        deduplicate_fastq.main(**vars(args))

    elif args.subcommand == 'digest_fastq':
        from capturec import digest_fastq
        del args.subcommand
        digest_fastq.main(**vars(args))










if __name__ == '__main__':
    run()
