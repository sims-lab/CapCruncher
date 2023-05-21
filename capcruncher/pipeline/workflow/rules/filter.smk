def validate_custom_filtering():

    custom_filter_stages = config["analysis"].get("custom_filtering", "")
    if not custom_filter_stages:
        cf = ""
    elif not os.path.exists(custom_filter_stages):
        cf = ""
    else:
        cf = f"--custom-filtering {custom_filter_stages}"

    return cf


rule filter_alignments:
    input:
        bam=rules.align_bowtie2.output.bam,
        annotations=rules.annotate.output.annotated,
    output:
        filtered_slices=temp(
            "capcruncher_output/alignment_filtering/initial/{sample}/{sample}_part{part}_{combined}.slices.parquet"
        ),
        stats_read="capcruncher_output/statistics/filtering/data/{sample}_part{part}_{combined}.read.stats.csv",
        stats_slice="capcruncher_output/statistics/filtering/data/{sample}_part{part}_{combined}.slice.stats.csv",
    params:
        analysis_method=config["analysis"]["method"],
        sample_name=lambda wildcards, output: wildcards.sample,
        output_prefix=lambda wildcards, output: output.filtered_slices.replace(
            ".slices.parquet", ""
        ),
        stats_prefix=lambda wildcards, output: output.stats_read.replace(
            ".read.stats.csv", ""
        ),
        read_type=lambda wildcards, output: wildcards.combined,
        custom_filtering=validate_custom_filtering(),
    log:
        "capcruncher_output/statistics/filtering/logs/{sample}_part{part}_{combined}.log",
    shell:
        """
        capcruncher \
        alignments \
        filter \
        {params.analysis_method} \
        -b {input.bam} \
        -a {input.annotations} \
        -o {params.output_prefix} \
        --stats-prefix {params.stats_prefix} \
        --sample-name {params.sample_name} \
        --read-type {params.read_type} \
        --no-cis-and-trans-stats \
        --no-fragments \
        {params.custom_filtering} > {log} 2>&1
        """


rule count_identified_viewpoints:
    input:
        slices=lambda wildcards: expand(
            "capcruncher_output/alignment_filtering/initial/{sample}/{sample}_part{part}_{combined}.slices.parquet",
            sample=[
                wildcards.sample,
            ],
            part=get_fastq_partition_numbers_for_sample(wildcards),
            combined=[
                "flashed",
                "pe",
            ],
        ),
    output:
        stats="capcruncher_output/statistics/identified_viewpoints/data/{sample}.identified_viewpoints.stats.csv",
    params:
        slices_dir = lambda wc: "capcruncher_output/alignment_filtering/initial/{sample}/",
    script:
        "../scripts/count_identified_viewpoints.py"


rule split_flashed_and_pe_datasets:
    input:
        slices_flashed=lambda wildcards: expand(
            "capcruncher_output/alignment_filtering/initial/{sample}/{sample}_part{part}_{combined}.slices.parquet",
            sample=[
                wildcards.sample,
            ],
            part=get_fastq_partition_numbers_for_sample(wildcards),
            combined=[
                "flashed",
            ],
        ),
        slices_pe=lambda wildcards: expand(
            "capcruncher_output/alignment_filtering/initial/{sample}/{sample}_part{part}_{combined}.slices.parquet",
            sample=[
                wildcards.sample,
            ],
            part=get_fastq_partition_numbers_for_sample(wildcards),
            combined=[
                "pe",
            ],
        ),
    output:
        slices_flashed=directory(
            "capcruncher_output/alignment_filtering/repartitioned/{sample}/flashed/"
        ),
        slices_pe=directory(
            "capcruncher_output/alignment_filtering/repartitioned/{sample}/pe/"
        ),
    shell:
        """
        mkdir -p {output.slices_flashed}
        mkdir -p {output.slices_pe}
        mv {input.slices_flashed} {output.slices_flashed}
        mv {input.slices_pe} {output.slices_pe}
        """


rule remove_duplicate_coordinates:
    input:
        slices_directory="capcruncher_output/alignment_filtering/repartitioned/{sample}/{combined}/",
    output:
        slices=directory(
            temp(
                "capcruncher_output/alignment_filtering/deduplicated/{sample}/{combined}"
            )
        ),
        stats_read="capcruncher_output/statistics/deduplication_by_coordinate/data/{sample}_{combined}.read.stats.csv",
    params:
        sample_name=lambda wildcards, output: wildcards.sample,
        stats_prefix=lambda wildcards, output: output.stats_read.replace(
            ".read.stats.csv", ""
        ),
        read_type=lambda wildcards, output: wildcards.combined,
    log:
        "logs/remove_duplicate_coordinates/{sample}_{combined}.log",
    script:
        "../scripts/remove_duplicate_coordinates.py"


rule combine_flashed_and_pe_post_deduplication:
    input:
        slices=expand(
            "capcruncher_output/alignment_filtering/deduplicated/{{sample}}/{combined}",
            combined=["flashed", "pe"],
        ),
    output:
        slices=directory("capcruncher_output/{sample}/{sample}.parquet"),
    shell:
        """
        mkdir -p {output.slices}
        mv {input.slices[0]}/*.parquet {output.slices}
        mv {input.slices[1]}/*.parquet {output.slices}
        """


rule cis_and_trans_stats:
    input:
        slices="capcruncher_output/{sample}/{sample}.parquet",
    output:
        stats="capcruncher_output/statistics/cis_and_trans_reporters/data/{sample}.reporter.stats.csv",
    params:
        sample_name=lambda wildcards, output: wildcards.sample,
        analysis_method=config["analysis"]["method"],
    log:
        "logs/cis_and_trans_stats/{sample}.log",
    shell:
        """
        capcruncher \
        utilities \
        cis-and-trans-stats \
        {input.slices} \
        --sample-name {params.sample_name} \
        """
