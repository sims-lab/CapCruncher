import json
import pyranges as pr


def validate_blacklist(blacklist):
    """Validate blacklist file."""

    blacklist_ok = True

    if blacklist is None:
        blacklist_ok = False
    elif not os.path.exists(blacklist):
        blacklist_ok = False

    return blacklist_ok


def configure_annotation_parameters():
    """Load defaults from annotation_defaults.json and overwrite with the current files"""
    defaults_file = workflow.source_path("../data/annotation_defaults.json")
    parameters = json.load(open(defaults_file))

    # Overwrite defaults with current options
    parameters["viewpoints"]["fn"] = config["analysis"]["viewpoints"]
    parameters["viewpoints"]["fraction"] = config["analysis_optional"].get(
        "minimum_viewpoint_overlap", parameters["viewpoints"]["fraction"]
    )
    parameters["viewpoints_count"]["fn"] = config["analysis"]["viewpoints"]
    parameters["viewpoints_count"]["fraction"] = config["analysis_optional"].get(
        "minimum_viewpoint_overlap", parameters["viewpoints"]["fraction"]
    )

    # Check if blacklist is valid
    if validate_blacklist(config["analysis"].get("blacklist")):
        parameters["blacklist"] = config["analysis"]["blacklist"]
    else:
        del parameters["blacklist"]

    return parameters


def format_annotation_parameters():
    """Format annotation parameters for use in the shell script."""

    parameters = configure_annotation_parameters()

    flags = {
        "name": "-n",
        "fn": "-b",
        "action": "-a",
        "fraction": "-f",
        "dtype": "-t",
    }

    annotation_args = []
    for annotation, options in parameters.items():
        for option, value in options.items():
            if value is not None:
                annotation_args.append(f"{flags[option]} {value}")

    return " ".join(annotation_args)


def format_priority_chromosome_list():
    """Format priority chromosome list for use in the shell script."""

    priority_chroms = config["analysis_optional"].get("priority_chromosomes", "")

    if not priority_chroms or priority_chroms == "None":
        chromosomes = None
    elif "," in priority_chroms:
        chromosomes = priority_chroms
    elif "viewpoints" in priority_chroms:
        pr_viewpoints = pr.read_bed(config["analysis"]["viewpoints"])
        chromosomes = ",".join(pr_viewpoints.Chromosome)

    return f"--priority-chroms {chromosomes}" if chromosomes else ""


def get_threads_for_annotation(annotation_files_and_params):
    return annotation_files_and_params.count("-b")


rule exclusions:
    input:
        viewpoints=config["analysis"]["viewpoints"],
    output:
        exclusions="capcruncher_output/annotate/exclude.bed",
    params:
        genome=config["genome"]["chrom_sizes"],
        exclusion_zone=config["analysis"]["reporter_exclusion_zone"],
    shell:
        """
        bedtools slop -i {input.viewpoints} -g {params.genome} -b {params.exclusion_zone} |
        bedtools subtract -a - -b  {input.viewpoints} > {output.exclusions}
        """


rule annotate:
    input:
        bam=rules.align_bowtie2.output.bam,
        exclusions="capcruncher_output/annotate/exclude.bed",
        viewpoints=config["analysis"]["viewpoints"],
    output:
        annotated="capcruncher_output/annotate/{sample}/{sample}_part{part}_{combined}.parquet",
    params:
        annotation_files_and_params=format_annotation_parameters(),
        priority_chromosomes=format_priority_chromosome_list(),
        prioritize_cis_slices="--prioritize-cis-slices"
        if config["analysis_optional"].get("prioritize_cis_slices", "")
        else "",
    threads: 6
    resources:
        mem_mb=5000,
    log:
        "capcruncher_output/logs/capcruncher_output/annotate/{sample}/{sample}_part{part}_{combined}.log",
    shell:
        """
        capcruncher \
        alignments \
        annotate \
        {input.bam} \
        -o \
        {output.annotated} \
        {params.annotation_files_and_params} \
        {params.priority_chromosomes} \
        {params.prioritize_cis_slices} \
        -p {threads} > {log} 2>&1
        """
