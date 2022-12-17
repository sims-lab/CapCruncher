# def aggregate(wc):
#     files = expand(
#         "capcruncher_filter/{sample}/slices/{sample}_part{part}_{combined}.parquet",
#         sample=wc.sample,
#         part=get_fastq_partition_numbers_for_sample(wc),
#         combined=["flashed", "pe"],
#     )

#     assert files
#     return files

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
            "capcruncher_filter/01_filtered/{sample}/{sample}_part{part}_{combined}.slices.parquet"
        ),
        stats_read="capcruncher_statistics/filtering/data/{sample}_part{part}_{combined}.read.stats.csv",
        stats_slice="capcruncher_statistics/filtering/data/{sample}_part{part}_{combined}.slice.stats.csv",
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
        "capcruncher_statistics/filtering/logs/{sample}_part{part}_{combined}.log",
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


rule split_flashed_and_pe_datasets:
    input:
        slices_flashed=lambda wildcards: expand(
            "capcruncher_filter/01_filtered/{sample}/{sample}_part{part}_{combined}.slices.parquet",
            sample=[
                wildcards.sample,
            ],
            part=get_fastq_partition_numbers_for_sample(wildcards),
            combined=[
                "flashed",
            ],
        ),
        slices_pe=lambda wildcards: expand(
            "capcruncher_filter/01_filtered/{sample}/{sample}_part{part}_{combined}.slices.parquet",
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
            "capcruncher_filter/02_split_flashed_and_pe/{sample}/flashed/"
        ),
        slices_pe=directory("capcruncher_filter/02_split_flashed_and_pe/{sample}/pe/"),
    shell:
        """
        mkdir -p {output.slices_flashed}
        mkdir -p {output.slices_pe}
        mv {input.slices_flashed} {output.slices_flashed}
        mv {input.slices_pe} {output.slices_pe}
        """


rule remove_duplicate_coordinates_flashed:
    input:
        slices_directory=rules.split_flashed_and_pe_datasets.output.slices_flashed,
    output:
        slices=temp(
            directory(
                "capcruncher_filter/03_remove_duplicate_coordinates/{sample}/flashed"
            )
        ),
    log:
        "logs/remove_duplicate_coordinates_flashed/{sample}.log"
    shell:
        """
        capcruncher \
        reporters \
        deduplicate \
        {input.slices_directory} \
        -o {output.slices} \
        --read-type flashed \
        > {log} 2>&1
        """


rule remove_duplicate_coordinates_pe:
    input:
        slices_directory=rules.split_flashed_and_pe_datasets.output.slices_pe,
    output:
        slices=temp(directory("capcruncher_filter/03_remove_duplicate_coordinates/{sample}/pe")),
    log:
        "logs/remove_duplicate_coordinates_pe/{sample}.log"
    shell:
        """
        capcruncher \
        reporters \
        deduplicate \
        {input.slices_directory} \
        -o {output.slices} \
        --read-type pe \
        > {log} 2>&1
        """


rule combine_flashed_and_pe_post_deduplication:
    input:
        slices_flashed=rules.remove_duplicate_coordinates_flashed.output.slices,
        slices_pe=rules.remove_duplicate_coordinates_pe.output.slices,
    output:
        slices=directory("capcruncher_filter/04_reporters/{sample}/"),
    shell:
        """
        mkdir -p {output.slices}
        mv {input.slices_flashed}/*.parquet {output.slices}
        mv {input.slices_pe}/*.parquet {output.slices}
        """
