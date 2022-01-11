import click
from numpy import multiply


@click.group()
def cli():
    """Contains methods for interaction counting, storing, bedgraph generation, comparisons."""


@cli.command()
@click.argument("union_bedgraph")
@click.option(
    "-n",
    "--capture_name",
    help="Name of capture probe, must be present in viewpoint file.",
    required=True,
)
@click.option(
    "-c",
    "--capture_viewpoints",
    help="Path to capture viewpoints bed file",
    required=True,
)
@click.option(
    "-o",
    "--output_prefix",
    help="Output prefix for pairwise statistical comparisons",
    default="out",
)
@click.option(
    "--design_matrix",
    help="Path tsv file containing sample annotations (N_SAMPLES * N_INFO_COLUMNS)",
    default=None,
)
@click.option(
    "--grouping_col", help="Column to use for grouping replicates", default="condition"
)
@click.option(
    "--threshold_count",
    help="Minimum count required to be considered for analysis",
    default=20,
    type=click.FLOAT,
)
@click.option(
    "--threshold_q",
    help="Upper threshold of q-value required for output.",
    default=0.05,
    type=click.FLOAT,
)
@click.option(
    "--threshold_mean",
    help="Minimum mean count required for output.",
    default=0,
    type=click.FLOAT,
)
def differential(*args, **kwargs):
    """
    Identifies differential interactions between conditions.

    Parses a union bedgraph containg reporter counts from at least two conditions with
    two or more replicates for a single capture probe and outputs differential interaction
    results. Following filtering to ensure that the number of interactions is above the required
    threshold (--threshold_count), diffxpy is used to run a wald test after
    fitting a negative binomial model to the interaction counts.The options to filter
    results can be filtered by a minimum mean value (threshold_mean) and/or
    maximum q-value (threshold-q) are also provided.

    Notes:

     Currently both the capture viewpoints and the name of the probe being analysed must
     be provided in order to correctly extract cis interactions.

     If a N_SAMPLE * METADATA design matrix has not been supplied, the script
     assumes that the standard replicate naming structure has been followed
     i.e. SAMPLE_CONDITION_REPLICATE_(1|2).fastq.gz.
    """

    from capcruncher.cli.reporters_differential import differential

    differential(*args, **kwargs)


@cli.command()
@click.argument("uri")
@click.option(
    "-n",
    "--viewpoint_names",
    help="Viewpoint to extract and convert to bedgraph, if not provided will transform all.",
    multiple=True,
)
@click.option("-o", "--output_prefix", help="Output prefix for bedgraphs")
@click.option(
    "--normalisation",
    help="Method to use interaction normalisation",
    default="raw",
    type=click.Choice(["raw", "n_cis", "region"]),
)
@click.option(
    "--normalisation-regions",
    help="Regions to use for interaction normalisation. The --normalisation method MUST be 'region'",
    default=None,
    type=click.STRING,
)
@click.option(
    "--binsize",
    help="Binsize to use for converting bedgraph to evenly sized genomic bins",
    default=0,
)
@click.option("--gzip", help="Compress output using gzip", default=False, is_flag=True)
@click.option(
    "--scale_factor",
    help="Scale factor to use for bedgraph normalisation",
    default=1e6,
    type=click.INT,
)
@click.option(
    "--sparse/--dense",
    help="Produce bedgraph containing just positive bins (sparse) or all bins (dense)",
    default=True,
)
@click.option(
    "-f",
    "--format",
    help="Output file format",
    type=click.Choice(["bedgraph", "bigwig"], case_sensitive=False),
    default="bedgraph",
)
def pileup(*args, **kwargs):
    """
    Extracts reporters from a capture experiment and generates a bedgraph file.

    Identifies reporters for a single probe (if a probe name is supplied) or all capture
    probes present in a capture experiment HDF5 file.

    The bedgraph generated can be normalised by the number of cis interactions for
    inter experiment comparisons and/or extract pilups binned into even genomic windows.
    """

    from capcruncher.cli.reporters_pileup import pileup

    pileup(*args, **kwargs)


@cli.command()
@click.argument("reporters")
@click.option("-o", "--output", help="Name of output file", default="CC_cooler.hdf5")
@click.option(
    "--remove_exclusions",
    default=False,
    help="Prevents analysis of fragments marked as proximity exclusions",
    is_flag=True,
)
@click.option(
    "--remove_capture",
    default=False,
    help="Prevents analysis of capture fragment interactions",
    is_flag=True,
)
@click.option(
    "--subsample",
    default=0,
    help="Subsamples reporters before analysis of interactions",
    type=float,
)
@click.option(
    "--low-memory",
    is_flag=True,
    default=False,
    help="Will perform counting in batches specifed by the chunksize to save memory (less accurate)",
)
@click.option(
    "--chunksize",
    default=int(2e6),
    help="Number of records to process at once",
)
@click.option(
    "-f",
    "--fragment-map",
    help="Path to digested genome bed file",
)
@click.option(
    "-v",
    "--viewpoint-path",
    help="Path to viewpoints file",
)
@click.option(
    "--cooler-output",
    "output_as_cooler",
    help="Output counts in cooler format",
    is_flag=True,
)
@click.option(
    "-p",
    "--n-cores",
    default=1,
    help="Number of cores to use for counting.",
    type=int,
)
@click.option(
    "-t",
    "--file-type",
    help="File format for input",
    default="auto",
    type=click.Choice(["auto", "tsv", "hdf5", "parquet"], case_sensitive=False),
)
def count(*args, **kwargs):
    """
    Determines the number of captured restriction fragment interactions genome wide.

    Parses a reporter slices tsv and counts the number of unique restriction fragment
    interaction combinations that occur within each fragment.

    Options to ignore unwanted counts e.g. excluded regions or capture fragments are provided.
    In addition the number of reporter fragments can be subsampled if required.
    """
    
    if kwargs.get("output_as_cooler"):
        if not kwargs.get("fragment_map"):
            raise ValueError("Restriction fragment map must be provided for cooler output")
        elif not kwargs.get("viewpoint_path"):
            raise ValueError("Viewpoint path must be provided for cooler output")


    from capcruncher.cli.reporters_count import count

    count(*args, **kwargs)


@cli.group()
def store():
    """
    Store reporter counts.

    These commands store and manipulate reporter restriction fragment interaction
    counts as cooler formated groups in HDF5 files.

    See subcommands for details.

    """


@store.command(name="fragments")
@click.argument("counts", required=True)
@click.option(
    "-f",
    "--fragment-map",
    help="Path to digested genome bed file",
    required=True,
)
@click.option(
    "-v",
    "--viewpoint-path",
    help="Path to viewpoints file",
    required=True,
)
@click.option(
    "-n",
    "--viewpoint-name",
    help="Name of viewpoint to store",
    default="",
)
@click.option(
    "-g",
    "--genome",
    help="Name of genome",
)
@click.option(
    "--suffix",
    help="Suffix to append after the capture name for the output file",
)
@click.option(
    "-o",
    "--output",
    help="Name of output file. (Cooler formatted hdf5 file)",
    default="out.hdf5",
)
def store_fragments(*args, **kwargs):
    """
    Stores restriction fragment interaction combinations at the restriction fragment level.

    Parses reporter restriction fragment interaction counts produced by
    "capcruncher reporters count" and gerates a cooler formatted group in an HDF5 File.
    See `https://cooler.readthedocs.io/en/latest/` for further details.
    """
    from capcruncher.cli.reporters_store import fragments

    fragments(*args, **kwargs)


@store.command(name="bins")
@click.argument("cooler_path", required=True)
@click.option(
    "-b",
    "--binsizes",
    help="Binsizes to use for windowing",
    default=(5000,),
    multiple=True,
    type=click.INT,
)
@click.option(
    "--normalise",
    is_flag=True,
    help="Enables normalisation of interaction counts during windowing",
)
@click.option(
    "--overlap_fraction",
    help="Minimum overlap between genomic bins and restriction fragments for overlap",
    default=0.5,
)
@click.option(
    "-p",
    "--n_cores",
    help="Number of cores used for binning",
    default=4,
    type=click.INT,
)
@click.option(
    "--scale_factor",
    help="Scaling factor used for normalisation",
    default=1e6,
    type=click.INT,
)
@click.option(
    "--conversion_tables",
    help="Pickle file containing pre-computed fragment -> bin conversions.",
    default=None,
)
@click.option(
    "-o",
    "--output",
    help="Name of output file. (Cooler formatted hdf5 file)",
    default="out.hdf5",
)
def store_bins(*args, **kwargs):
    """
    Convert a cooler group containing restriction fragments to constant genomic windows

    Parses a cooler group and aggregates restriction fragment interaction counts into
    genomic bins of a specified size. If the normalise option is selected,
    columns containing normalised counts are added to the pixels table of the output
    """
    from capcruncher.cli.reporters_store import bins

    bins(*args, **kwargs)


@store.command(name="merge")
@click.argument("coolers", required=True, nargs=-1)
@click.option("-o", "--output", help="Output file name")
def store_merge(*args, **kwargs):
    """
    Merges capcruncher cooler files together.

    Produces a unified cooler with both restriction fragment and genomic bins whilst
    reducing the storage space required by hard linking the "bins" tables to prevent duplication.
    """
    from capcruncher.cli.reporters_store import merge

    merge(*args, **kwargs)


@cli.group()
def compare():

    """
    Compare bedgraphs and CapCruncher cooler files.

    These commands allow for specific viewpoints to be extracted from cooler files
    and perform user defined groupby aggregations.

    See subcommands for details.

    """


@compare.command(name="concat")
@click.argument("infiles", required=True, nargs=-1)
@click.option(
    "-f",
    "--format",
    help="Input file format",
    type=click.Choice(["auto", "bedgraph", "cooler"]),
    default="cooler",
)
@click.option("-o", "--output", help="Output file name", default="union.tsv")
@click.option("-v", "--viewpoint", help="Viewpoint to extract")
@click.option("-r", "--resolution", help="Resolution to extract")
@click.option(
    "--region", help="Limit to specific coordinates in the format chrom:start-end"
)
@click.option(
    "--normalisation",
    help="Method to use interaction normalisation",
    default="raw",
    type=click.Choice(["raw", "n_cis", "region"]),
)
@click.option(
    "--normalisation-regions",
    help="Regions to use for interaction normalisation. The --normalisation method MUST be 'region'",
    default=None,
    type=click.STRING,
)
@click.option(
    "--scale_factor",
    help="Scale factor to use for bedgraph normalisation",
    default=1e6,
    type=click.INT,
)
@click.option(
    "-p", "--n_cores", help="Number of cores to use for extracting bedgraphs", default=1
)
def bedgraphs_concat(*args, **kwargs):

    from capcruncher.cli.reporters_compare import concat

    concat(*args, **kwargs)


@compare.command(name="summarise")
@click.argument("infile", required=True)
@click.option("-o", "--output-prefix", help="Output file prefix")
@click.option(
    "-f", "--output-format", type=click.Choice(["bedgraph", "tsv"]), default="bedgraph"
)
@click.option(
    "-m",
    "--summary-methods",
    help="Summary methods to use for aggregation. Can be any method in numpy or scipy.stats",
    multiple=True,
)
@click.option("-n", "--group-names", help="Group names for aggregation", multiple=True)
@click.option(
    "-c",
    "--group-columns",
    help="Column names/numbers (0 indexed, the first column after the end coordinate counts as 0) for aggregation.",
    multiple=True,
)
@click.option(
    "--subtraction", is_flag=True, help="Perform subtration between aggregated groups"
)
@click.option("--suffix", help="Add a suffix before the file extension")
def bedgraphs_summarise(*args, **kwargs):

    from capcruncher.cli.reporters_compare import summarise

    summarise(*args, **kwargs)
