def get_count_files(wc):
    counts = []
    counts.append(
        f"capcruncher_output/pileups/counts_by_restriction_fragment/{wc.sample}.hdf5"
    )

    if PERFORM_BINNING:
        counts.append(
            f"capcruncher_output/pileups/counts_by_genomic_bin/{wc.sample}.hdf5"
        )

    return counts


def get_normalisation_from_config(wc):
    regions = config["normalisation"]["regions"]

    if not regions is None or isinstance(regions, str):
        if os.path.exists(regions):
            return f"--normalisation region --normalisation-regions {regions}"
    return "--normalisation n_cis"


if CAPCRUNCHER_TOOLS:

    rule count:
        input:
            slices=rules.combine_flashed_and_pe_post_deduplication.output.slices,
            restriction_fragment_map=rules.digest_genome.output.bed,
            viewpoints=VIEWPOINTS,
        output:
            temp(
                "capcruncher_output/pileups/counts_by_restriction_fragment/{sample}.hdf5"
            ),
        log:
            "capcruncher_output/logs/counts/{sample}.log",
        threads: 8
        resources:
            mem_mb=lambda wc, attempt: 3000 * 2**attempt,
        params:
            outdir="capcruncher_output/pileups/counts_by_restriction_fragment",
        shell:
            """
            mkdir -p {params.outdir} && \
            capcruncher-tools \
            count \
            {input.slices} \
            -o {output} \
            -f {input.restriction_fragment_map} \
            -v {input.viewpoints} \
            -p {threads} \
            > {log} 2>&1
            """

else:

    rule count:
        input:
            slices=rules.combine_flashed_and_pe_post_deduplication.output.slices,
            restriction_fragment_map=rules.digest_genome.output.bed,
            viewpoints=VIEWPOINTS,
        output:
            temp(
                "capcruncher_output/pileups/counts_by_restriction_fragment/{sample}.hdf5"
            ),
        log:
            "capcruncher_output/logs/counts/{sample}.log",
        threads: 8
        resources:
            mem_mb=lambda wc, attempt: 3000 * 2**attempt,
        params:
            outdir="capcruncher_output/pileups/counts_by_restriction_fragment",
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
            --cooler-output \
            -p {threads} \
            > {log} 2>&1
            """


rule bin_counts:
    input:
        "capcruncher_output/pileups/counts_by_restriction_fragment/{sample}.hdf5",
    output:
        temp("capcruncher_output/pileups/counts_by_genomic_bin/{sample}.hdf5"),
    params:
        bin_size=[f"-b {b}" for b in BIN_SIZES],
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
        > {log} 2>&1
        """


rule merge_counts:
    input:
        get_count_files,
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
        bedgraph="capcruncher_output/interim/pileups/bedgraphs/{sample}/raw/{sample}_{viewpoint}.bedgraph",
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
        bedgraph="capcruncher_output/interim/pileups/bedgraphs/{sample}/norm/{sample}_{viewpoint}.bedgraph",
    log:
        "capcruncher_output/logs/bedgraph_norm/{sample}_{viewpoint}.log",
    params:
        output_prefix=lambda wc, output: pathlib.Path(output.bedgraph).parent
        / f"{wc.sample}",
        viewpoint=lambda wc, output: wc.viewpoint,
        normalisation=get_normalisation_from_config,
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
    log:
        "capcruncher_output/logs/bedgraph_to_bigwig/{sample}_{norm}_{viewpoint}.log",
    params:
        chrom_sizes=config["genome"]["chrom_sizes"],
    shell:
        """
        sort -k1,1 -k2,2n {input.bedgraph} > {input.bedgraph}.sorted
        bedGraphToBigWig {input.bedgraph}.sorted {params.chrom_sizes} {output.bigwig}
        rm {input.bedgraph}.sorted
        """
