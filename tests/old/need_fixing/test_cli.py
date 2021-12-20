# import glob
# import os

# import click
# import pysam
# import pytest
# import ujson
# import xopen
# from click.testing import CliRunner
# from capcruncher.cli import cli
# import pandas as pd


# # Pre-run setup
# dir_test = os.path.realpath(os.path.dirname(__file__))
# dir_package = os.path.dirname(dir_test)
# dir_data = os.path.join(dir_package, "data")


# @pytest.fixture(scope="session", autouse=True)
# def cleanup():

#     yield

#     # Remove previous test files
#     for fn in glob.glob(os.path.join(dir_test, "test", "*")):
#         os.unlink(fn)

#     for fn in glob.glob(os.path.join(dir_test, "stats", "*")):
#         os.unlink(fn)


# def get_fastq_n_records(fq):
#     return sum(1 for l in pysam.FastxFile(fq))

# def test_cli_runs():
#     """Test checks that the cli is functional and the help option works"""

#     runner = CliRunner()
#     result = runner.invoke(cli, ["--help"])
#     assert result.exit_code == 0


# def test_genome_digest():

#     test_fa = os.path.join(dir_data, "test", "genome_digest", "test.fa")
#     test_output = os.path.join(dir_test, "test", "test_digest_genome.bed")
#     test_output_stats = os.path.join(dir_test, "stats", "test_digest_genome.bed")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "genome",
#             "digest",
#             test_fa,
#             "-r",
#             "dpnii",
#             "-o",
#             test_output,
#             "-l",
#             test_output_stats,
#             "--sort",
#         ],
#     )

#     assert result.exit_code == 0

#     with open(test_output) as r:
#         test_bed = [l.strip().split("\t") for l in r]

#     # Should be just the one cutsite at the moment
#     assert len(test_bed) == 2


# def test_fastq_deduplicate_parse():

#     # Test parsing
#     output_parsed = os.path.join(dir_test, "test", "fq_parsed_test.json")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "fastq",
#             "deduplicate",
#             "parse",
#             os.path.join(dir_data, "test",  "fastq_deduplication","duplicated_1.fastq.gz"),
#             os.path.join(dir_data, "test",  "fastq_deduplication","duplicated_2.fastq.gz"),
#             "-o",
#             output_parsed,
#         ],
#     )

#     # Check that the script exits successfully
#     assert result.exit_code == 0

#     with open(output_parsed) as f:
#         result_test = ujson.load(f)

#     with open(os.path.join(dir_data, "test",  "fastq_deduplication","fq_parsed_result.json")) as f:
#         result_correct = ujson.load(f)

#     # Check that all keys and values match between the saved file and the test
#     assert result_test == result_correct


# def test_fastq_deduplicate_identification():
#     output_duplicates_test = os.path.join(dir_test, "test", "fq_duplicates_test.json")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "fastq",
#             "deduplicate",
#             "identify",
#             os.path.join(dir_data, "test",  "fastq_deduplication","fq_parsed_result.json"),
#             "-o",
#             output_duplicates_test,
#         ],
#     )

#     assert result.exit_code == 0

#     with open(output_duplicates_test) as f:
#         result_test = ujson.load(f)

#     with open(os.path.join(dir_data, "test",  "fastq_deduplication","fq_duplicates_results.json")) as f:
#         result_correct = ujson.load(f)

#     # Checks that the same number of duplicates are identified, some randomness in identification
#     assert len(result_test) == len(result_correct)


# def test_fastq_deduplicate_removal():

#     # Test removal
#     output_removal_prefix = os.path.join(dir_test, "test", "fq_dedup_test")
#     output_removal_test = output_removal_prefix + "_1.fastq"
#     output_removal_test_stats = os.path.join(dir_test, "stats", "deduplication")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "fastq",
#             "deduplicate",
#             "remove",
#             "-d",
#             os.path.join(dir_data, "test",  "fastq_deduplication","fq_duplicates_results.json"),
#             "-o",
#             output_removal_prefix,
#             os.path.join(dir_data, "test",  "fastq_deduplication","duplicated_1.fastq.gz"),
#             os.path.join(dir_data, "test",  "fastq_deduplication","duplicated_2.fastq.gz"),
#             "--stats-prefix",
#             output_removal_test_stats,
#             "--sample-name",
#             "test",
#         ],
#     )

#     assert result.exit_code == 0
#     assert os.path.exists(output_removal_test_stats + ".deduplication.csv")

#     with open(output_removal_test) as r:
#         result_test = r.readlines()

#     with open(os.path.join(dir_data, "test",  "fastq_deduplication","fq_dedup_1.fastq")) as r:
#         result_correct = r.readlines()

#     with open(os.path.join(dir_data, "test",  "fastq_deduplication","fq_duplicates_results.json")) as r:
#         duplicates = ujson.load(r)

#     fq_unfilt_n_entries = get_fastq_n_records(os.path.join(dir_data, "test",  "fastq_deduplication","duplicated_1.fastq.gz"))
#     fq_dd_n_entries = get_fastq_n_records(output_removal_test)

#     # Checks the number of expected duplicates are removed
#     assert len(duplicates) == fq_unfilt_n_entries - fq_dd_n_entries

#     # Checks the test file matches the expected file
#     assert len(result_test) == len(result_correct)


# def test_fastq_digest():

#     fq1 = os.path.join(dir_data, "test", "fastq_digestion","digest_1.fastq.gz")
#     fq2 = os.path.join(dir_data, "test", "fastq_digestion" ,"digest_2.fastq.gz")
#     test_output_flashed = os.path.join(
#         dir_test, "test", "test_fastq_digest_flashed.fastq"
#     )
#     test_output_pe = os.path.join(dir_test, "test", "test_fastq_digest_pe.fastq")
#     test_output_stats = os.path.join(dir_test, "stats", "digestion")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "fastq",
#             "digest",
#             "-m",
#             "pe",
#             "-r",
#             "dpnii",
#             "-o",
#             test_output_pe,
#             "--stats-prefix",
#             test_output_stats,
#             "--sample-name",
#             "test",
#             fq1,
#             fq2,
#         ],
#     )

#     assert result.exit_code == 0
#     assert os.path.exists(test_output_stats + ".digestion.read.summary.csv")

#     os.unlink(test_output_stats + ".digestion.read.summary.csv")
#     result = runner.invoke(
#         cli,
#         [
#             "fastq",
#             "digest",
#             "-m",
#             "flashed",
#             "-r",
#             "dpnii",
#             "-o",
#             test_output_flashed,
#             "--stats-prefix",
#             test_output_stats,
#             "--sample-name",
#             "test",
#             fq1,
#         ],
#     )

#     assert result.exit_code == 0

#     # check stats are produced
#     assert os.path.exists(test_output_stats + ".digestion.read.summary.csv")


# def test_alignments_annotate():

#     test_output = os.path.join(dir_test, "test", "test_annotate.hdf5")
#     test_bed = os.path.join(dir_data, "test", "alignment_annotation", "test_slices.bed")
#     test_capture = os.path.join(dir_data, "test", "alignment_annotation","test_capture.bed")
#     test_exclusion = os.path.join(dir_data, "test", "alignment_annotation","test_exclusions.bed")
#     test_rf = os.path.join(dir_data, "test", "alignment_annotation","test_rf.bed")
#     test_bad_bed = os.path.join(dir_data, "test", "alignment_annotation", "bad_bed.bed")

#     # Test with bad exclusions file
#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "alignments",
#             "annotate",
#             test_bed,
#             "-b",
#             test_bad_bed,
#             "-n",
#             "will_fail",
#             "-a",
#             "count"
#         ],
#     )

#     assert result.exit_code == 1

#     # Test with ignoring bad file
#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "alignments",
#             "annotate",
#             test_bed,
#             "-b",
#             test_capture,
#             "-b",
#             test_capture,
#             "-b",
#             test_exclusion,
#             "-b",
#             test_rf,
#             "-a",
#             "get",
#             "-a",
#             "count",
#             "-a",
#             "get",
#             "-a",
#             "get",
#             "-n",
#             "capture",
#             "-n",
#             "capture_count",
#             "-n",
#             "exclusion",
#             "-n",
#             "rf",
#             "-o",
#             test_output,
#             "--invalid_bed_action",
#             "ignore",
#         ],
#     )

#     assert result.exit_code == 0

#     # Test with bad path
#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "alignments",
#             "annotate",
#             "XXXXXXXXXXXX",
#         ],
#     )

#     assert result.exit_code != 0


# def test_alignments_filter():

#     bam = os.path.join(dir_data, "test", "Slc25A37-test_1_part0.pe.bam")
#     annotations = os.path.join(dir_data, "test", "Slc25A37-test_1_part0.pe.annotations.tsv")
#     output_prefix = os.path.join(dir_test, "test", "test_filtering_cli")
#     stats_prefix = os.path.join(dir_test, "stats", "test_filtering_cli")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "alignments",
#             "filter",
#             "capture",
#             "-b",
#             bam,
#             "-a",
#             annotations,
#             "-o",
#             output_prefix,
#             "--stats-prefix",
#             stats_prefix,
#             "--fragments",
#             "--read-stats",
#             "--slice-stats",
#             "--cis-and-trans-stats",
#         ],
#     )

#     assert result.exit_code == 0

#     df_slices = pd.read_csv(f"{output_prefix}.Slc25A37.slices.tsv", sep="\t")
#     assert df_slices.shape[0] == 1062


# def test_reporter_deduplicate_identify():
    
#     reporters_duplicated = os.path.join(dir_data, "test", "reporters_duplicated.tsv")
#     reporters_duplicated_ids = os.path.join(dir_test, "test", "reporters_deduplicated.json")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "alignments",
#             "deduplicate",
#             "identify",
#             reporters_duplicated,
#             '--read-type',
#             'flashed',
#             "-o",
#             reporters_duplicated_ids
#         ],
#     )

#     assert result.exit_code == 0
#     with xopen.xopen(reporters_duplicated_ids, "r") as r:
#         ids_duplicated = {int(x) for x in ujson.load(r)}
#     assert len(ids_duplicated) == 2

#     result = runner.invoke(
#         cli,
#         [
#             "alignments",
#             "deduplicate",
#             "identify",
#             reporters_duplicated,
#             '--read-type',
#             'pe',
#             "-o",
#             reporters_duplicated_ids
#         ],
#     )
    
#     assert result.exit_code == 0
#     with xopen.xopen(reporters_duplicated_ids, "r") as r:
#         ids_duplicated = {int(x) for x in ujson.load(r)}
#     assert len(ids_duplicated) == 2






# def test_reporter_count():

#     reporters = os.path.join(dir_data, "test", "Slc25A37_reporters.tsv.gz")
#     output = os.path.join(dir_test, "test", "test_reporter_counts.tsv")
#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "count",
#             reporters,
#             "-o",
#             output,
#         ],
#     )

#     assert result.exit_code == 0

#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "count",
#             reporters,
#             "-o",
#             output,
#             "--low-memory",
#         ],
#     )

#     assert result.exit_code == 0

#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "count",
#             reporters,
#             "-o",
#             output,
#             "--subsample",
#             "0.8",
#             "--remove_exclusions",
#             "--remove_capture",
#         ],
#     )

#     assert result.exit_code == 0


# def test_reporter_storage():

#     counts = os.path.join(dir_data, "test", "Slc25A37_reporter_counts.tsv.gz")
#     bins = os.path.join(dir_data, "test", "genome.digest.bed.gz")
#     viewpoints = os.path.join(dir_data, "test", "mm9_capture_oligos.bed")
#     output_prefix = os.path.join(dir_test, "test/cli_cooler")
#     output = os.path.join(dir_test, "test/cli_cooler.Slc25A37.fragments.hdf5")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "store",
#             "fragments",
#             counts,
#             "-f",
#             bins,
#             "-g",
#             "mm9",
#             "-n",
#             "Slc25A37",
#             "-o",
#             output_prefix,
#             "-v",
#             viewpoints,
#             "--suffix",
#             "fragments",
#         ],
#     )

#     infile = output
#     output_prefix = os.path.join(dir_test, "test", "cli_cooler_binned.hdf5")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "store",
#             "bins",
#             output,
#             "-b",
#             "2500",
#             "-o",
#             output_prefix,
#             "--normalise",
#             "-p",
#             "4",
#         ],
#     )

#     assert result.exit_code == 0

#     runner = CliRunner()
#     clrs = glob.glob(os.path.join(dir_test, "test", "cli*.hdf5"))


#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "store",
#             "merge",
#             *clrs,
#             '-o',
#             os.path.join(dir_test, "test", "cli_cooler_merged.hdf5"),
#         ],
#     )

#     assert result.exit_code == 0

# def test_reporters_pileup():
    
#     runner = CliRunner()
#     clr_path = os.path.join(dir_data, "test", "test.Slc25A37.hdf5")
#     bdg_prefix = os.path.join(dir_test, "test", "test.Slc25A37")
#     bdg_path = f"{bdg_prefix}..bedgraph"


#     result = runner.invoke(
#         cli,
#         [
#             "reporters",
#             "pileup",
#             clr_path,
#             '-o',
#             bdg_prefix,
#         ],
#     )

#     assert result.exit_code == 0
#     assert os.path.exists(bdg_path)
 






# def test_plot_make_templates():

#     try:
#         import coolbox
#         template = os.path.join(dir_test, "test", "test_pileup_template.yml")
#         runner = CliRunner()
        
#         result = runner.invoke(
#             cli,
#             [
#                 "plot",
#                 "make-template",
#                 *glob.glob(os.path.join(dir_data, "test", "test_bigwigs", "*.bigWig")),
#                 os.path.join(dir_data, "test", "mm9_chr14_genes.bed"),
#                 "-b",
#                 "5000",
#                 "-d",
#                 os.path.join(dir_data, "test", "design_matrix.tsv"),
#                 "-o",
#                 template.replace(".yml", "")
#             ],
#         )
        
#         assert result.exit_code == 0
#         assert os.path.exists(template)
    
#     except ImportError:
#         pass
    


# def test_plot_make_plots():

#     try:
#         import coolbox
#         template = os.path.join(dir_test, "test", "test_pileup_template.yml")
#         plot = os.path.join(dir_test, "test", "test_pileup.svg")

#         runner = CliRunner()
        
#         result = runner.invoke(
#             cli,
#             [
#                 "plot",
#                 "make-plot",
#                 "-c",
#                 template,
#                 "-r",
#                 "chr14:69878554-69933221",
#                 "-o",
#                 plot,
#                 "--x-axis",
#             ],
#         )
        
#         assert result.exit_code == 0
#         assert os.path.exists(plot)

#     except ImportError:
#         pass



# def test_gtf_to_bed():


#     gtf = os.path.join(dir_data, "test", "mm9_chr14.gtf")
#     bed = os.path.join(dir_test, "test", "mm9_chr14.bed")

#     runner = CliRunner()
#     result = runner.invoke(
#         cli,
#         [
#             "utilities",
#             "gtf-to-bed12",
#             gtf,
#             "-o",
#             bed,
#         ],
#     )

#     assert result.exit_code == 0
#     assert os.path.exists(bed)
