import os
import sys
import argparse

# Make sure the script can find capturec scripts
SCRIPT_PATH = os.path.abspath(__file__)
SCRIPT_DIR = os.path.dirname(SCRIPT_PATH)
PACKAGE_DIR = os.path.dirname(SCRIPT_DIR)
sys.path.append(PACKAGE_DIR)


def prepare_parser():
    description = "Collection of scripts to analyse Capture-C, Tri-C and Tiled-C data."
    parser = argparse.ArgumentParser(
        prog="ccanalyser",
        description=description,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    return parser


def add_category_args(subcommand):

    categories = sorted(["utils", "validate", "ccanalysis", "stats"])
    return {cat: subcommand.add_parser(cat) for cat in categories}


def add_agg_stats_args(subcommand):

    parser = subcommand.add_parser("agg_stats", help="Combine run statistics")

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


def add_aggregate_tsv_args(subcommand):
    parser = subcommand.add_parser("agg_tsv", help="Combine run statistics")
    subparser = parser.add_subparsers(dest="method")

    parser_join = subparser.add_parser("join")
    parser_join.add_argument("-i", "--input_files", nargs="+", required=True)
    parser_join.add_argument("--index", default="parent_read", required=True)
    parser_join.add_argument("--header", default=False, action="store_true")
    parser_join.add_argument("-o", "--output", default="joined.tsv.gz")
    parser_join.add_argument("-p", "--n_processes", default=8, type=int)

    parser_concatenate = subparser.add_parser("concatenate")
    parser_concatenate.add_argument("-i", "--input_files", nargs="+", required=True)
    parser_concatenate.add_argument("--header", default=False, action="store_true")
    parser_concatenate.add_argument("-o", "--output", default="concatenated.tsv.gz")

    parser_aggregate = subparser.add_parser("aggregate")
    parser_aggregate.add_argument("-i", "--input_files", required=True)
    parser_aggregate.add_argument("--index", default=None)
    parser_aggregate.add_argument("--header", default=False, action="store_true")
    parser_aggregate.add_argument("-o", "--output", default="aggregated.tsv.gz")
    parser_aggregate.add_argument("-g", "--groupby_columns", nargs="+")
    parser_aggregate.add_argument("--aggregate_method", nargs="+")
    parser_aggregate.add_argument("--aggregate_columns", nargs="+")


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

    # TODO: Sort out duplicate options
    parser.add_argument("--duplicates", default="remove", choices=["remove",])


def add_ccanalyser_args(subcommand):
    parser = subcommand.add_parser(
        "ccanalysis", help="Analysis of Capture-C processed BAM file"
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


def add_slices_to_bedgraph_args(subcommand):

    parser = subcommand.add_parser(
        "slices_to_bdg", help="Converts ccanalyser output to bedgraph"
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


def add_count_interactions_args(subcommand):
    parser = subcommand.add_parser(
        "count_interactions", help="Counts restriction fragment combinations"
    )
    parser.add_argument("-f", "--slices")
    parser.add_argument("-o", "--outfile", default="out.tsv.gz")
    parser.add_argument(
        "--only_cis",
        default=False,
        action="store_true",
        help="Only count cis interactions",
    )
    parser.add_argument(
        "--remove_exclusions",
        default=False,
        action="store_true",
        help="Remove proximity exclusions",
    )
    parser.add_argument("--subsample", default=0, type=int, help="Subsample fragments")


def add_deduplicate_fastq_args(subcommand):

    parser = subcommand.add_parser(
        "deduplicate_fastq", help="Removes PCR duplicates from fastq file"
    )

    subparser = parser.add_subparsers(dest="mode")
    parser_fd = subparser.add_parser("find_duplicates")
    parser_fd.add_argument('input_files', nargs='+')
    parser_fd.add_argument('-d', '--deduplicated_ids', required=True)
    parser_fd.add_argument(
        "-c",
        "--compression_level",
        help="Level of compression (1-9 with 9 being the highest)",
        type=int,
        default=5,
    )
    parser_fd.add_argument(
        "--read_buffer",
        help="defines the number of reads processed before writing to file",
        default=100000,
        type=int,
    )

    parser_md = subparser.add_parser("merge_ids")
    parser_md.add_argument('input_files', nargs='+')
    parser_md.add_argument('-o', '--output_files', default='merged.pkl')


    parser_rd = subparser.add_parser("remove_duplicates")
    parser_rd.add_argument('input_files', nargs='+')
    parser_rd.add_argument('-d', '--deduplicated_ids', required=True, nargs='+')
    parser_rd.add_argument('-o', '--output_files', required=True, nargs='+')
    parser_rd.add_argument(
        "-c",
        "--compression_level",
        help="Level of compression (1-9 with 9 being the highest)",
        type=int,
        default=5,
    )
    parser_rd.add_argument(
        "--read_buffer",
        help="defines the number of reads processed before writing to file",
        default=100000,
        type=int,
    )




def add_digest_fastq_args(subcommand):

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "-o", "--output_file", help="output file name", default="digested.fastq.gz",
    )

    parent_parser.add_argument(
        "-r", "--restriction_enzyme", help="Name or sequence of restriction enzyme"
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

    parser.add_argument(
        "-r", "--recognition_site", help="Name of restriction enzyme or its sequence"
    )
    parser.add_argument(
        "-l", "--logfile", help="filename for logfile", default="test.log"
    )


def add_plot_matrix_args(subcommand):
    parser = subcommand.add_parser("plot_interactions", help="Plots interaction matrix")
    parser.add_argument(
        "cooler_restriction_fragments",
        help="HDF5 file in Cooler format containing counts of restriction fragment combination counts",
    )
    parser.add_argument(
        "cooler_binned",
        help="HDF5 file in Cooler format containing binned counts of restriction fragment combination counts",
    )
    parser.add_argument("-c", "--coordinates", required=True, type=str)
    parser.add_argument("-m", "--method", choices=["tri", "tiled"])
    parser.add_argument(
        "-n",
        "--normalisation",
        default=None,
        choices=["infer", "raw", "ice", "scaling_factor"],
    )
    parser.add_argument("--scaling_factor", default=1e5, type=int)
    parser.add_argument("--binning_method", default="overlap")
    parser.add_argument("--overlap_fraction", default=0.5, type=float)
    parser.add_argument("--filter_low_counts", default=0.02, type=float)
    parser.add_argument("--filter_high_counts", default=0, type=float)
    parser.add_argument("-o", "--output_prefix", default="")
    parser.add_argument(
        "-f", "--output_format", default="png", choices=["png", "svg", "jpeg"]
    )
    parser.add_argument("--cmap", help="Colour map to use", default="viridis")


def add_remove_duplicate_slices_args(subcommand):
    parser = subcommand.add_parser("rmdup_slices", help="Removes duplicate slices")
    subparser = parser.add_subparsers(dest="mode")

    parser_fragments = subparser.add_parser("fragments")
    parser_fragments.add_argument("-i", "--input_files", nargs="+", required=True)
    parser_fragments.add_argument("-f", "--deduplicated_fragments", required=True)
    parser_fragments.add_argument(
        "--shuffle",
        help="shuffles the input files to randomise the deduplication",
        action="store_true",
    )
    parser_fragments.add_argument("-p", "--n_cores", default=8, type=int)
    parser_fragments.add_argument("-m", "--max_memory", default="64GB", type=str)

    parser_slices = subparser.add_parser("slices")
    parser_slices.add_argument("-i", "--input_files", required=True)
    parser_slices.add_argument("-f", "--deduplicated_fragments", required=True)
    parser_slices.add_argument("-o", "--output", default="deduplicated.tsv.gz")
    parser_slices.add_argument("-p", "--n_cores", default=8, type=int)
    parser_slices.add_argument("-m", "--max_memory", default="64GB", type=str)


def add_split_fastq_args(subcommand):

    parser = subcommand.add_parser(
        "split_fastq", help="Splits fastq file into smaller chunks"
    )
    
    parser.add_argument("input_files", help="Fasq file(s) to parse", nargs='+')

    parser.add_argument(
        "-n", '--n_reads', help="Number of reads per output file", default=1000000, type=int
    )

    parser.add_argument(
        "-c",
        "--compression_level",
        help="Level of gzip compression (1-9 with 9 being the most compressed/slowest)",
        default=6,
        type=int,
    )
    parser.add_argument("-o", "--output_prefix", help="output prefix", default="split")


def add_store_interactions_args(subcommand):
    parser = subcommand.add_parser(
        "store_interactions",
        help="Stores restriction fragment interactions as .cool files",
    )
    parser.add_argument("-c", "--restriction_fragment_counts", required=True)
    parser.add_argument("-m", "--restriction_fragment_map", required=True)
    parser.add_argument("-g", "--genome", default="mm9", required=True, type=str)
    parser.add_argument("-b", "--binsizes", default=1000, nargs="+", type=int)
    parser.add_argument("--bin_method", default="overlap")
    parser.add_argument("-f", "--overlap_fraction", default=0.51, type=float)
    parser.add_argument("-o", "--output_prefix", default="counts.hdf5")


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
            add_split_fastq_args,
            add_aggregate_tsv_args,
        ],
        "ccanalysis": [
            add_ccanalyser_args,
            add_slices_to_bedgraph_args,
            add_annotate_slices_args,
            add_count_interactions_args,
            add_store_interactions_args,
            add_remove_duplicate_slices_args,
            add_plot_matrix_args
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
            from ccanalyser.utils import digest_fastq

            digest_fastq.main(**vars(args))

        elif subcommand == "digest_genome":
            from ccanalyser.utils import digest_genome

            digest_genome.main(**vars(args))

        elif subcommand == "agg_tsv":
            from ccanalyser.utils import aggregate_tsv

            aggregate_tsv.main(**vars(args))

        elif subcommand == "split_fastq":
            from ccanalyser.utils import split_fastq

            split_fastq.main(**vars(args))

    elif category == "ccanalysis":

        if subcommand == "ccanalysis":
            from ccanalyser.ccanalysis import ccanalysis

            ccanalysis.main(**vars(args))

        elif subcommand == "slices_to_bdg":
            from ccanalyser.ccanalysis import slices_to_bedgraph

            slices_to_bedgraph.main(**vars(args))

        elif subcommand == "annotate_slices":
            from ccanalyser.ccanalysis import annotate_slices

            annotate_slices.main(**vars(args))

        elif subcommand == "count_interactions":
            from ccanalyser.ccanalysis import count_interactions

            count_interactions.main(**vars(args))

        elif subcommand == "store_interactions":
            from ccanalyser.ccanalysis import store_interactions

            store_interactions.main(**vars(args))
        
        elif subcommand == 'rmdup_slices':
            from ccanalyser.ccanalysis import remove_duplicate_slices

            remove_duplicate_slices.main(**vars(args))
        
        elif subcommand == 'plot_interactions':
            from ccanalyser.ccanalysis import plot_interactions

            plot_interactions.main(**vars(args))


    elif category == "validate":
        if subcommand == "test_parameters":
            from ccanalyser.validate import test_run_params

            test_run_params.main()

        elif subcommand == "validate_annotations":
            from ccanalyser.validate import validate_annotations

            validate_annotations.main(**vars(args))

    elif category == "stats":
        if subcommand == "agg_stats":
            from ccanalyser.stats import aggregate_statistics

            aggregate_statistics.main(**vars(args))


if __name__ == "__main__":
    main()
