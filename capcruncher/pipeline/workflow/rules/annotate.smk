import json
import pyranges as pr
import capcruncher.pipeline.utils


rule exclusions:
    input:
        viewpoints=config["analysis"]["viewpoints"],
    output:
        exclusions="capcruncher_output/interim/annotate/exclude.bed",
    params:
        genome=config["genome"]["chrom_sizes"],
        exclusion_zone=config["analysis"]["reporter_exclusion_zone"],
    shell:
        """
        bedtools slop -i {input.viewpoints} -g {params.genome} -b {params.exclusion_zone} |
        bedtools subtract -a - -b  {input.viewpoints} > {output.exclusions}
        """


rule check_n_bins_per_viewpoint:
    input:
        bins=rules.digest_genome.output.bed,
        viewpoints=config["analysis"]["viewpoints"],
    output:
        sentinel="capcruncher_output/resources/validation/check_n_bins_per_viewpoint.sentinel",
        n_bins_per_viewpoint="capcruncher_output/resources/validation/n_bins_per_viewpoint.tsv",
    params:
        ignore_multiple_bins_per_viewpoint=IGNORE_MULTIPLE_FRAGMENTS_PER_VIEWPOINT
        if ASSAY in ["capture", "tri"]
        else True,
    log:
        "capcruncher_output/logs/validation/check_n_bins_per_viewpoint.log",
    script:
        "../scripts/validation_check_n_bins_per_viewpoint.py"


rule annotate:
    input:
        bam=rules.align_bowtie2.output.bam,
        exclusions="capcruncher_output/interim/annotate/exclude.bed",
        viewpoints=config["analysis"]["viewpoints"],
        single_bin_per_viewpoint=rules.check_n_bins_per_viewpoint.output.sentinel,
    output:
        annotated=temp(
            "capcruncher_output/interim/annotate/{sample}/{sample}_part{part}_{combined}.parquet"
        ),
    params:
        annotation_files_and_params=capcruncher.pipeline.utils.format_annotation_parameters(
            workflow, config
        ),
        priority_chromosomes=capcruncher.pipeline.utils.format_priority_chromosome_list(
            config
        ),
        prioritize_cis_slices="--prioritize-cis-slices"
        if config["analysis_optional"].get("prioritize_cis_slices", "")
        else "",
    threads: 6
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * 2**attempt,
    log:
        "capcruncher_output/logs/annotate/{sample}/{sample}_part{part}_{combined}.log",
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
