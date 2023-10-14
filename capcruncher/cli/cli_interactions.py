import click


@click.group()
def cli():
    """Contains methods for interaction counting, storing, bedgraph generation, comparisons."""


@cli.command()
@click.argument("slices", required=True)
@click.option(
    "-o",
    "--output",
    help="Output prefix for directory of deduplicated slices",
    default="deduplicated_slices/",
)
@click.option(
    "--statistics",
    help="Output prefix for stats file(s)",
    default="",
)
@click.option(
    "--sample-name", help="Name of sample e.g. DOX_treated_1", default="sample"
)
@click.option(
    "--read-type",
    help="Type of read",
    default="flashed",
    type=click.Choice(["flashed", "pe"], case_sensitive=False),
)
def deduplicate(*args, **kwargs):
    """
    Identifies and removes duplicated aligned fragments.

    PCR duplicates are very commonly present in Capture-C/Tri-C/Tiled-C data and must be removed
    for accurate analysis. Unlike fastq deduplicate, this command removes fragments with identical
    genomic coordinates.

    Non-combined (pe) and combined (flashed) reads are treated slightly differently due to the increased
    confidence that the ligation junction has been captured for the flashed reads.

    """

    from capcruncher.cli.interactions_deduplicate import deduplicate

    deduplicate(*args, **kwargs)


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
    "--scale-factor",
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

    from capcruncher.cli.interactions_pileup import pileup

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
    "-p",
    "--n-cores",
    default=1,
    help="Number of cores to use for counting.",
    type=int,
)
@click.option(
    "--assay", type=click.Choice(["capture", "tri", "tiled"]), default="capture"
)
def count(*args, **kwargs):
    """
    Determines the number of captured restriction fragment interactions genome wide.

    Counts the number of interactions between each restriction fragment and all other
    restriction fragments in the fragment.

    The output is a cooler formatted HDF5 file containing a single group containing
    the interactions between restriction fragments.

    See `https://cooler.readthedocs.io/en/latest/` for further details.

    """

    from capcruncher.cli.interactions_count import count

    count(*args, **kwargs)


@cli.command(name="counts-to-cooler")
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
    from capcruncher.cli.interactions_store import fragments

    fragments(*args, **kwargs)


@cli.command(name="fragments-to-bins")
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
    "--scale-factor",
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
@click.option(
    "--assay", type=click.Choice(["capture", "tri", "tiled"]), default="capture"
)
def store_bins(*args, **kwargs):
    """
    Convert a cooler group containing restriction fragments to constant genomic windows

    Parses a cooler group and aggregates restriction fragment interaction counts into
    genomic bins of a specified size. If the normalise option is selected,
    columns containing normalised counts are added to the pixels table of the output
    """
    from capcruncher.cli.interactions_store import bins

    bins(*args, **kwargs)


@cli.command(name="merge")
@click.argument("coolers", required=True, nargs=-1)
@click.option("-o", "--output", help="Output file name")
def store_merge(*args, **kwargs):
    """
    Merges capcruncher HDF5 files together.

    Produces a unified cooler with both restriction fragment and genomic bins whilst
    reducing the storage space required by hard linking the "bins" tables to prevent duplication.
    """
    from capcruncher.api.storage import merge_coolers

    merge_coolers(*args, **kwargs)


@cli.group()
def compare():

    r"""Compare bedgraphs and CapCruncher cooler files.

    These commands allow for specific viewpoints to be extracted from CapCruncher HDF5 files and perform:

        1. User defined groupby aggregations.

        2. Comparisons between conditions.

        3. Identification of differential interactions between conditions.

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

    from capcruncher.cli.interactions_compare import concat

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

    from capcruncher.cli.interactions_compare import summarise

    summarise(*args, **kwargs)


@compare.command(name="differential")
@click.argument("interaction_files", required=True, nargs=-1)
@click.option(
    "-o", "--output-prefix", help="Output file prefix", default="differential"
)
@click.option("-v", "--viewpoint", help="Viewpoint to extract", required=True)
@click.option("-d", "--design-matrix", help="Design matrix file", required=True)
@click.option("-c", "--contrast", help="Contrast to test", default="condition")
@click.option(
    "-r",
    "--regions-of-interest",
    help="Regions of interest to test for differential interactions",
    default=None,
)
@click.option(
    "--viewpoint-distance",
    help="Distance from viewpoint to test for differential interactions",
    default=None,
    type=click.INT,
)
@click.option(
    "--threshold-count",
    help="Minimum number of interactions to test for differential interactions",
    default=20,
)
@click.option(
    "--threshold-q",
    help="Minimum q-value to test for differential interactions",
    default=0.05,
)
def bedgraphs_differential(*args, **kwargs):
    """Perform differential testing on CapCruncher HDF5 files.

    This command performs differential testing on CapCruncher HDF5 files. It requires a design matrix
    and a contrast to test. The design matrix should be a tab separated file with the first column
    containing the sample names and the remaining columns containing the conditions. The contrast
    should specify the name of the column in the design matrix to test. The output is a tab separated bedgraph.
    """
    from capcruncher.cli.interactions_differential import differential

    differential(*args, **kwargs)
