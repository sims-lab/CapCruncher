import os
import sys
import argparse
import imp

# Make sure the script can find capturec scripts
SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(PACKAGE_DIR)


def prepare_parser():
    description = "Collection of scripts to analyse NG-Capture-C data."
    parser = argparse.ArgumentParser(
        prog="capturec",
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    return parser


def add_category_args(subcommand):

    categories = sorted(["utils", "validate", "ccanalysis", "stats"])
    return {cat: subcommand.add_parser(cat) for cat in categories}


def add_agg_stats_args(subcommand):

    parser = subcommand.add_parser("aggregate_stats", help="Combine run statistics")

    parser.add_argument(
        "--deduplication_stats",
        nargs="+",
        help="Deduplication stats paths",
        required=True,
    )
    parser.add_argument(
        "--digestion_stats", nargs="+", help="Digestion stats paths", required=True
    )
    parser.add_argument("--ccanalyser_stats", nargs="+", required=True)
    parser.add_argument("--reporter_stats", nargs="+", required=True)
    parser.add_argument("--output_dir", default="run_stats")


def add_annotate_slices_args(subcommand):

    parser = subcommand.add_parser(
        "annotate_slices", help="Intersects bed files and generates dataframe"
    )

    parser.add_argument("-a", "--bed1")
    parser.add_argument("-b", "--bed2", nargs="+")
    parser.add_argument("--actions", nargs="+", choices=["get", "count"])
    parser.add_argument("-c", "--colnames", nargs="+")
    parser.add_argument(
        "-f", "--overlap_fractions", nargs="*", default=1e-9, type=float
    )
    parser.add_argument("-o", "--outfile", default="out.tsv.gz")

    #TODO: Sort out duplicate options
    parser.add_argument('--duplicates', default='remove', choices=['remove',])


def add_ccanalyser_args(subcommand):
    parser = subcommand.add_parser(
        "ccanalyser", help="Analysis of Capture-C processed BAM file"
    )
    parser.add_argument("-i", "--input_bam", help="BAM file to parse", required=True)
    parser.add_argument(
        "-a",
        "--annotations",
        help="Tab-delimited text file containing annotation for each read in bam file",
        required=True,
    )
    parser.add_argument(
        "--output_prefix", help="Output file prefix", default="ccanalyser_out"
    )
    parser.add_argument(
        "--stats_output", help="stats files output prefix", default="stats"
    )

    parser.add_argument(
        "-m",
        "--method",
        help="Determines the mode to run",
        default="capture",
        choices=["capture", "tri", "tiled"],
    )


def add_convert_tsv_to_bedgraph_args(subcommand):

    parser = subcommand.add_parser(
        "ccanalyser_to_bedgraph", help="Converts ccanalyser output to bedgraph"
    )
    parser.add_argument("-i", "--tsv_input", help="Reporter tsv file")
    parser.add_argument(
        "-b",
        "--bed",
        help="""Bed file to intersect with reporters
                                      e.g. RE fragments bed file.""",
    )
    parser.add_argument(
        "-o", "--output", help="Output file name", default="out.bedgraph"
    )


def add_count_rf_combs_args(subcommand):
    parser = subcommand.add_parser(
        "count_rf_combs", help="Counts restriction fragment combinations"
    )
    parser.add_argument("-f", "--slices")
    parser.add_argument("-o", "--outfile", default="out.tsv.gz")
    parser.add_argument('--only_cis', default=False, action='store_true', help='Only count cis interactions')


def add_deduplicate_fastq_args(subcommand):

    parser = subcommand.add_parser(
        "deduplicate_fastq", help="Removes PCR duplicates from fastq file"
    )

    parser.add_argument(
        "-1", "--fq1", help="fastq file to parse containing read 1", required=True
    )
    parser.add_argument(
        "-2", "--fq2", help="fastq file to parse containing read 2", required=True
    )
    parser.add_argument(
        "--out1", help="fastq file to parse containing read 1", default="out1.fq.gz"
    )
    parser.add_argument(
        "--out2", help="fastq file to parse containing read 2", default="out2.fq.gz"
    )
    parser.add_argument(
        "--stats_file", help="name of deduplication statistics file", default=sys.stdout
    )
    parser.add_argument(
        "-c",
        "--compression_level",
        help="Level of compression (1-9 with 9 being the highest)",
        type=int,
        default=5,
    )
    parser.add_argument(
        "--read_buffer",
        help="defines the number of reads processed before writing to file",
        default=10000,
        type=int,
    )


def add_digest_fastq_args(subcommand):

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "-o", "--output_file", help="output file name", default="digested.fastq.gz",
    )

    enzyme_group = parent_parser.add_mutually_exclusive_group(required=True)
    enzyme_group.add_argument(
        "-r", "--restriction_enzyme", help="Name of restriction enzyme"
    )
    enzyme_group.add_argument(
        "-s", "--cut_sequence", help="Sequence of restriction site"
    )

    parent_parser.add_argument(
        "-m",
        "--minimum_slice_length",
        help="Shortest length for a slice to be output",
        default=20,
        type=int,
    )
    parent_parser.add_argument(
        "--stats_file", help="stats_file_prefix", default="stats.log"
    )
    parent_parser.add_argument(
        "-c",
        "--compression_level",
        help="Level of gzip compression (1-9 with 9 being the most compressed/slowest)",
        default=6,
        type=int,
    )
    parent_parser.add_argument(
        "--keep_cutsite",
        help="Determines if cutsite is stripped from the start of each slice",
        action="store_true",
        default=False,
    )
    parent_parser.add_argument(
        "--buffer",
        help="Number of reads to process before writing output",
        default=10000,
        type=int,
    )

    parent_parser.add_argument(
        "-p",
        "--n_digestion_processes",
        help="Number of digestion processes to spawn",
        default=1,
        type=int,
    )

    parser = subcommand.add_parser(
        "digest_fastq", help="Performs in silico digestion of fastq files",
    )
    subparsers = parser.add_subparsers(
        help="Run in either flashed or unflashed", dest="subcommand"
    )

    parser_flashed = subparsers.add_parser(
        "flashed", help="For flashed reads", parents=[parent_parser,]
    )
    parser_flashed.add_argument(
        "-i", "--input_fastq", help="fastq file to parse", required=True
    )

    parser_unflashed = subparsers.add_parser(
        "unflashed", help="For unflashed reads", parents=[parent_parser,]
    )
    parser_unflashed.add_argument(
        "-1", "--fq1", help="fastq file containing read 1 to parse", required=True
    )
    parser_unflashed.add_argument(
        "-2", "--fq2", help="fastq file containing read 2 to parse", required=True
    )


def add_digest_genome_args(subcommand):

    parser = subcommand.add_parser(
        "digest_genome", help="Identifies all restriction fragments in a given genome"
    )
    parser.add_argument(
        "-i", "--input_fasta", help="fasta file to parse", required=True
    )
    parser.add_argument(
        "-o", "--output_file", help="output file name", default="digested.bed"
    )

    enzyme_group = parser.add_mutually_exclusive_group(required=True)
    enzyme_group.add_argument(
        "-r", "--restriction_enzyme", help="Name of restriction enzyme"
    )
    enzyme_group.add_argument(
        "-s", "--cut_sequence", help="Sequence of restriction site"
    )
    parser.add_argument(
        "-l", "--logfile", help="filename for logfile", default="test.log"
    )


def add_join_tsv_args(subcommand):
    parser = subcommand.add_parser(
        "join_tsv", help="Concatenates or joins multiple tsv files"
    )
    parser.add_argument(
        "-f",
        "--index_field",
        help="shared column name to join on",
        nargs="?",
        default=None,
    )
    parser.add_argument(
        "-o", "--output_file", help="output file name", default="joined.tsv"
    )
    parser.add_argument("-i", "--input_files", nargs="+", help="at least 2 input files")
    parser.add_argument(
        "-m",
        "--method",
        help="join method to use (join/concatenate)",
        choices=["join", "concatenate"],
        default="join",
    )


def add_split_fastq_args(subcommand):

    parser = subcommand.add_parser(
        "split_fastq", help="Splits fastq file into smaller chunks"
    )
    parser.add_argument("-i", "--input_fastq", help="BAM file to parse", required=True)
    parser.add_argument(
        "--chunksize", help="Number of reads per output file", default=1000000, type=int
    )
    parser.add_argument(
        "-c",
        "--compression_level",
        help="Level of gzip compression (1-9 with 9 being the most compressed/slowest)",
        default=6,
        type=int,
    )
    parser.add_argument("-n", "--output_prefix", help="output prefix", default="split")


def add_validation_args(subcommand):
    parser = subcommand.add_parser(
        "validate_annotations", help="Validates ccanalyser annotation input"
    )
    parser.add_argument("-i", "--input_file", help="input .tsv file")
    parser.add_argument(
        "-o", "--output_file", help="output file name", default="joined.tsv"
    )
    parser.add_argument(
        "--col_names", help="additional col names to expect", nargs="+", default=[]
    )
    parser.add_argument(
        "--dtypes", help="dtypes of additional columns", nargs="+", default=[]
    )


def add_test_parameters_args(subcommand):
    parser = subcommand.add_parser(
        "test_parameters",
        help="Confirms that parameters in capturec_pipeline.yml are valid run parameters",
    )


def main():

    parser = prepare_parser()
    subcommand = parser.add_subparsers(dest="category", required=True)
    category_parsers = add_category_args(subcommand)

    category_scripts = {
        "utils": [
            add_deduplicate_fastq_args,
            add_digest_genome_args,
            add_digest_fastq_args,
            add_join_tsv_args,
            add_split_fastq_args,
        ],
        "ccanalysis": [
            add_ccanalyser_args,
            add_convert_tsv_to_bedgraph_args,
            add_annotate_slices_args,
            add_count_rf_combs_args,
        ],
        "validate": [add_validation_args, add_test_parameters_args,],
        "stats": [add_agg_stats_args,],
    }

    for name, p in category_parsers.items():
        scripts_args = category_scripts[name]
        cmd = p.add_subparsers(dest="command", required=True)

        for script_arg in scripts_args:
            script_arg(cmd)

    args = parser.parse_args()
    category = args.category
    subcommand = args.command

    # Remove unexpected args from the namespace
    del args.category
    del args.command

    if category == "utils":

        if subcommand == "deduplicate_fastq":
            from ccanalyser.utils import deduplicate_fastq

            deduplicate_fastq.main(**vars(args))

        elif subcommand == "digest_fastq":
            from ccanalyser.utils import digest_fastq_dev

            digest_fastq_dev.main(**vars(args))

        elif subcommand == "digest_genome":
            from ccanalyser.utils import digest_genome

            digest_genome.main(**vars(args))

        elif subcommand == "join_tsv":
            from ccanalyser.utils import join_tsv

            join_tsv.main(**vars(args))

        elif subcommand == "split_fastq":
            from ccanalyser.utils import split_fastq

            split_fastq.main(**vars(args))

    elif category == "ccanalysis":

        if subcommand == "ccanalyser":
            from ccanalyser.ccanalysis import ccanalyser

            ccanalyser.main(**vars(args))

        elif subcommand == "ccanalyser_to_bedgraph":
            from ccanalyser.ccanalysis import convert_tsv_to_bedgraph

            convert_tsv_to_bedgraph.main(**vars(args))

        elif subcommand == "annotate_slices":
            from ccanalyser.ccanalysis import annotate_slices

            annotate_slices.main(**vars(args))

        elif subcommand == "count_rf_combs":
            from ccanalyser.ccanalysis import count_restriction_fragment_combinations

            count_restriction_fragment_combinations.main(**vars(args))

    elif category == "validate":
        if subcommand == "test_parameters":
            from ccanalyser.validate import test_run_params

            test_run_params.main()

        elif subcommand == "validate_annotations":
            from ccanalyser.validate import validate_annotations

            validate_annotations.main(**vars(args))

    elif category == "stats":
        if subcommand == "aggregate_stats":
            from ccanalyser.stats import aggregate_statistics

            aggregate_statistics.main(**vars(args))


if __name__ == "__main__":
    main()
