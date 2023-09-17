import capcruncher.pipeline.utils


def get_mem_mb(wildcards, threads, attempt=0):
    return threads * 3000 * 2 ** (attempt - 1)


def get_outdir(wildcards, output):
    return str(pathlib.Path(output[0]).parent)


rule count:
    input:
        slices=rules.combine_flashed_and_pe_post_deduplication.output.slices,
        restriction_fragment_map=rules.digest_genome.output.bed,
        viewpoints=VIEWPOINTS,
    output:
        temp(
            "capcruncher_output/interim/pileups/counts_by_restriction_fragment/{sample}.hdf5"
        ),
    log:
        "capcruncher_output/logs/counts/{sample}.log",
    threads: 8
    resources:
        mem_mb=get_mem_mb,
    params:
        outdir=get_outdir,
        assay=config["analysis"]["method"],
    shell:
        """
        mkdir -p {params.outdir} && \
        capcruncher \
        interactions \
        count \
        {input.slices} \
        -o {output} \
        -f {input.restriction_fragment_map} \
        -v {input.viewpoints} \
        -p {threads} \
        --assay {params.assay}
        > {log} 2>&1
        """


rule bin_counts:
    input:
        "capcruncher_output/interim/pileups/counts_by_restriction_fragment/{sample}.hdf5",
    output:
        temp("capcruncher_output/interim/pileups/counts_by_genomic_bin/{sample}.hdf5"),
    params:
        bin_size=[f"-b {b}" for b in BIN_SIZES],
        assay=config["analysis"]["method"],
    log:
        "capcruncher_output/logs/bin_counts/{sample}.log",
    threads: 4
    resources:
        mem_mb=lambda wc, attempt: 3000 * 2**attempt,
    shell:
        """
        capcruncher \
        interactions \
        fragments-to-bins \
        {input} \
        -o {output} \
        {params.bin_size} \
        -p {threads} \
        --assay {params.assay} \
        > {log} 2>&1
        """


rule merge_counts:
    input:
        lambda wc: capcruncher.pipeline.utils.get_count_files(wc, PERFORM_BINNING),
    output:
        "capcruncher_output/results/{sample}/{sample}.hdf5",
    log:
        "capcruncher_output/logs/merge_counts/{sample}.log",
    threads: 1
    shell:
        """
        capcruncher \
        interactions \
        merge \
        {input} \
        -o {output} \
        > {log} 2>&1
        """


rule bedgraph_raw:
    input:
        cooler=rules.merge_counts.output,
    output:
        bedgraph=temp(
            "capcruncher_output/interim/pileups/bedgraphs/{sample}/raw/{sample}_{viewpoint}.bedgraph"
        ),
    retries: 0
    log:
        "capcruncher_output/logs/bedgraph_raw/{sample}_{viewpoint}.log",
    params:
        output_prefix=lambda wc, output: pathlib.Path(output.bedgraph).parent
        / f"{wc.sample}",
        viewpoint=lambda wc, output: wc.viewpoint,
    shell:
        """
        capcruncher \
        interactions \
        pileup \
        {input.cooler} \
        -o {params.output_prefix} \
        -n {params.viewpoint} \
        --normalisation raw \
        > {log} 2>&1
        """


rule bedgraph_normalised:
    input:
        cooler=rules.merge_counts.output,
    output:
        bedgraph=temp(
            "capcruncher_output/interim/pileups/bedgraphs/{sample}/norm/{sample}_{viewpoint}.bedgraph"
        ),
    log:
        "capcruncher_output/logs/bedgraph_norm/{sample}_{viewpoint}.log",
    retries: 0
    params:
        output_prefix=lambda wc, output: pathlib.Path(output.bedgraph).parent
        / f"{wc.sample}",
        viewpoint=lambda wc, output: wc.viewpoint,
        normalisation=lambda wc: capcruncher.pipeline.utils.get_normalisation_from_config(
            wc, config
        ),
        scale_factor=config["normalisation"].get("scale_factor", int(1e6)),
    shell:
        """
        capcruncher \
        interactions \
        pileup \
        {input.cooler} \
        -o {params.output_prefix} \
        -n {params.viewpoint} \
        {params.normalisation} \
        --scale-factor {params.scale_factor} \
        > {log} 2>&1
        """


rule bedgraph_to_bigwig:
    input:
        bedgraph="capcruncher_output/interim/pileups/bedgraphs/{sample}/{norm}/{sample}_{viewpoint}.bedgraph",
    output:
        bigwig="capcruncher_output/results/{sample}/bigwigs/{norm}/{sample}_{viewpoint}.bigWig",
    retries: 0
    log:
        "capcruncher_output/logs/bedgraph_to_bigwig/{sample}_{norm}_{viewpoint}.log",
    params:
        chrom_sizes=config["genome"]["chrom_sizes"],
    shell:
        """
        sort -k1,1 -k2,2n {input.bedgraph} > {input.bedgraph}.sorted
        bedGraphToBigWig {input.bedgraph}.sorted {params.chrom_sizes} {output.bigwig} 2> {log}
        rm {input.bedgraph}.sorted
        """
