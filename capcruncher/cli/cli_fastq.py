import click
import pathlib
import ast


class OptionEatAll(click.Option):
    def __init__(self, *args, **kwargs):
        self.save_other_options = kwargs.pop("save_other_options", True)
        nargs = kwargs.pop("nargs", -1)
        assert nargs == -1, "nargs, if set, must be -1 not {}".format(nargs)
        super(OptionEatAll, self).__init__(*args, **kwargs)
        self._previous_parser_process = None
        self._eat_all_parser = None

    def add_to_parser(self, parser, ctx):
        def parser_process(value, state):
            # method to hook to the parser.process
            done = False
            value = [value]
            if self.save_other_options:
                # grab everything up to the next option
                while state.rargs and not done:
                    for prefix in self._eat_all_parser.prefixes:
                        if state.rargs[0].startswith(prefix):
                            done = True
                    if not done:
                        value.append(state.rargs.pop(0))
            else:
                # grab everything remaining
                value += state.rargs
                state.rargs[:] = []
            value = tuple(value)

            # call the actual process
            self._previous_parser_process(value, state)

        retval = super(OptionEatAll, self).add_to_parser(parser, ctx)
        for name in self.opts:
            our_parser = parser._long_opt.get(name) or parser._short_opt.get(name)
            if our_parser:
                self._eat_all_parser = our_parser
                self._previous_parser_process = our_parser.process
                our_parser.process = parser_process
                break
        return retval



@click.group()
def cli():
    """Contains methods for fastq splitting, deduplicating and digestion."""


@cli.command()
@click.argument("input_files", nargs=-1, required=True)
@click.option(
    "-m",
    "--method",
    help="Method to use for splitting",
    type=click.Choice(["python", "unix"]),
    default="unix",
)
@click.option(
    "-o",
    "--output_prefix",
    help="Output prefix for deduplicated fastq file(s)",
    default="split",
)
@click.option(
    "--compression_level",
    help="Level of compression for output files",
    default=5,
    type=click.INT,
)
@click.option(
    "-n",
    "--n_reads",
    help="Number of reads per fastq file",
    default=1e6,
    type=click.INT,
)
@click.option(
    "--gzip/--no-gzip", help="Determines if files are gziped or not", default=False
)
@click.option("-p", "--n_cores", default=1, type=click.INT)
@click.option(
    "-s",
    "--suffix",
    help="Suffix to add to output files (ignore {read_number}.fastq as this is added automatically)",
    default="",
)
def split(*args, **kwargs):
    """
    Splits fastq file(s) into equal chunks of n reads.

    """

    from capcruncher.cli.fastq_split import split

    split(*args, **kwargs)


def digest_fastq(
    fastqs: tuple,
    restriction_site: str,
    mode: str = "pe",
    output_file: pathlib.Path | str = "out.fastq.gz",
    minimum_slice_length: int = 18,
    statistics: pathlib.Path | str = "digest.json",
    sample_name: str = "sampleX",
    **kwargs,
):
    """Digest FASTQ files and write digestion statistics."""

    from capcruncher_tools.api import digest_fastq as digest_fastq_records
    from loguru import logger

    from capcruncher.utils import get_restriction_site

    logger.info("Digesting FASTQ files")

    if len(fastqs) > 1 and mode == "flashed":
        raise ValueError("Flashed mode can only be used with a single FASTQ file")

    stats = digest_fastq_records(
        fastqs=fastqs,
        restriction_site=get_restriction_site(restriction_site),
        output=output_file,
        read_type=mode.title(),
        sample_name=sample_name,
        minimum_slice_length=minimum_slice_length,
    )

    logger.info("Digestion complete. Generating statistics")
    with open(statistics, "w") as f:
        f.write(stats.model_dump_json())

    return stats


@cli.command()
@click.argument("fastqs", nargs=-1, required=True)
@click.option(
    "-r",
    "--restriction_enzyme",
    help="Restriction enzyme name or sequence to use for in silico digestion.",
    required=True,
)
@click.option(
    "-m",
    "--mode",
    help="Digestion mode. Combined (Flashed) or non-combined (PE) read pairs.",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
    required=True,
)
@click.option("-o", "--output_file", default="out.fastq.gz")
@click.option("--minimum_slice_length", default=18, type=click.INT)
@click.option("--statistics", help="Output path for stats file", default="stats")
@click.option(
    "--sample-name",
    help="Name of sample e.g. DOX_treated_1. Required for correct statistics.",
    default="sampleX",
)
def digest(*args, **kwargs):
    """
    Performs in silico digestion of one or a pair of fastq files.
    """
    kwargs["restriction_site"] = kwargs.pop("restriction_enzyme")
    digest_fastq(*args, **kwargs)


def deduplicate_fastq(
    fastq_1: list[str],
    fastq_2: list[str],
    output_prefix: str | pathlib.Path = "deduplicated_",
    statistics: str = "deduplication_statistics.json",
    sample_name: str = "sampleX",
    shuffle: bool = False,
    **kwargs,
):
    """Deduplicate paired FASTQ files and write deduplication statistics."""

    from capcruncher.api.statistics import FastqDeduplicationStatistics
    from capcruncher_tools.api import deduplicate_fastq as deduplicate_fastq_records
    from loguru import logger

    df_stats = deduplicate_fastq_records(
        fastq1=fastq_1,
        fastq2=fastq_2,
        output_prefix=output_prefix,
        sample_name=sample_name,
        shuffle=shuffle,
    )

    dedup_stats = FastqDeduplicationStatistics(
        sample=sample_name,
        total=df_stats.query("stat_type == 'read_pairs_total'")["stat"].values[0],
        duplicates=df_stats.query("stat_type == 'read_pairs_duplicated'")[
            "stat"
        ].values[0],
    )
    with open(statistics, "w") as f:
        f.write(dedup_stats.model_dump_json())

    logger.info("Printing deduplication statistics to stdout")
    df_vis = df_stats.copy()
    df_vis["stat_type"] = df_vis["stat_type"].str.replace("_", " ").str.title()
    df_vis = df_vis[["stat_type", "stat"]]
    df_vis.columns = ["Stat Type", "Number of Reads"]
    print(df_vis.to_string(index=False))


@cli.command()
@click.option(
    "-1", "--fastq1", help="Read 1 FASTQ files", required=True, cls=OptionEatAll
)
@click.option(
    "-2", "--fastq2", help="Read 2 FASTQ files", required=True, cls=OptionEatAll
)
@click.option(
    "-o",
    "--output-prefix",
    help="Output prefix for deduplicated FASTQ files",
    default="deduped",
)
@click.option(
    "--sample-name", help="Name of sample e.g. DOX_treated_1", default="sampleX"
)
@click.option(
    "-s", "--statistics", help="Statistics output file name", default="stats.csv"
)
@click.option(
    "--shuffle",
    help="Shuffle reads before deduplication",
    is_flag=True,
    default=False,
)
def deduplicate(*args, **kwargs):
    """
    Identifies PCR duplicate fragments from Fastq files.

    PCR duplicates are very commonly present in Capture-C/Tri-C/Tiled-C data and must be removed
    for accurate analysis. These commands attempt to identify and remove duplicate reads/fragments
    from fastq file(s) to speed up downstream analysis.

    """
    fq1 = [pathlib.Path(f) for f in ast.literal_eval(kwargs["fastq1"])]
    fq2 = [pathlib.Path(f) for f in ast.literal_eval(kwargs["fastq2"])]

    kwargs["fastq_1"] = fq1
    kwargs["fastq_2"] = fq2

    deduplicate_fastq(*args, **kwargs)
