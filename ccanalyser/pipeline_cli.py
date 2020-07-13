import argparse
import os
import sys
import importlib
import subprocess as sub

# Make sure the script can find capturec scripts
SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(PACKAGE_DIR)


def prepare_parser():
    description = 'Pipeline to automate processing of Capture-C data'
    parser = argparse.ArgumentParser(
        prog='capturec',
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    return parser


def add_get_yaml_args(subcommand):
    parser = subcommand.add_parser(
        'get_yaml', help='Gets an example copy of the parameters file.'
    )

def add_pipeline_args(subcommand):
    parser = subcommand.add_parser(
        'pipeline', help='Runs the analysis pipeline'
    )

def main():

    parser = prepare_parser()
    subcommand = parser.add_subparsers(dest='command')

    for func in [add_get_yaml_args, add_pipeline_args]:
        subparser = func(subcommand)

    args, raw_args = parser.parse_known_args()
    subcommand = args.command

    if subcommand == 'get_yaml':
        from ccanalyser.pipeline import get_yaml
    elif subcommand == 'pipeline':
        cmd = ['python',
               os.path.join(SCRIPT_DIR, 'pipeline', 'capturec_pipeline.py'),
               *raw_args]
        sub.run(cmd)

if __name__ == '__main__':
    main()
