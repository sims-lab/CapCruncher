def is_binning_required():

    bin_sizes = config["analysis"].get("bin_sizes", None)

    if (bin_sizes is not None) or bin_size != "None" or bin_size != 0:
        if isinstance(bin_sizes, int):
            return True
        elif isinstance(bin_sizes, list):
            if all(
                isinstance(bin_size, int) and bin_size > 0 for bin_size in bin_sizes
            ):
                return True


def get_count_files(wc):

    counts = []
    counts.append(f"capcruncher_pileup/01_counts_by_fragment/{wc.sample}.hdf5")

    binning_required = is_binning_required()

    if binning_required:
        counts.append(f"capcruncher_pileup/02_counts_by_genomic_bin/{wc.sample}.hdf5")

    return counts


def get_normalisation_from_config():
    regions = config["normalisation"]["regions"]

    if not regions is None or isinstance(regions, str):
        if os.path.exists(regions):
            return f"--normalisation region --normalisation-regions {regions}"
    return "--normalisation n_cis"


rule count:
    input:
        slices=rules.combine_flashed_and_pe_post_deduplication.output.slices,
        restriction_fragment_map=rules.digest_genome.output.bed,
    output:
        "capcruncher_pileup/01_counts_by_fragment/{sample}.hdf5",
    params:
        viewpoints=config["analysis"]["viewpoints"],
    log:
        "logs/counts/{sample}.log",
    threads: 4
    shell:
        """
        capcruncher \
        interactions \
        count \
        {input.slices} \
        -o {output} \
        -f {input.restriction_fragment_map} \
        -v {params.viewpoints} \
        --cooler-output \
        -p {threads} \
        > {log} 2>&1
        """


rule bin_counts:
    input:
        "capcruncher_pileup/01_counts_by_fragment/{sample}.hdf5",
    output:
        "capcruncher_pileup/02_counts_by_genomic_bin/{sample}.hdf5",
    params:
        bin_size=config["analysis"]["bin_sizes"],
    log:
        "logs/bin_counts/{sample}.log",
    threads: 4
    shell:
        """
        capcruncher \
        interactions \
        fragments-to-bins \
        {input} \
        -o {output} \
        -b {params.bin_size} \
        -p {threads} \
        > {log} 2>&1
        """


rule merge_counts:
    input:
        get_count_files,
    output:
        "capcruncher_pileup/03_counts/{sample}.hdf5",
    log:
        "logs/merge_counts/{sample}.log",
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
        bedgraph=expand(
            "capcruncher_pileup/04_bedgraph/{{sample}}/raw/{{sample}}.raw.{viewpoint}.bedgraph",
            viewpoint=VIEWPOINT_NAMES,
        ),
    log:
        "logs/bedgraph_raw/{sample}.log",
    params:
        output_prefix=lambda wc, output: f"capcruncher_pileup/04_bedgraph/{wc.sample}/raw/{wc.sample}.raw",
    shell:
        """
        capcruncher \
        interactions \
        pileup \
        {input.cooler} \
        -o {params.output_prefix} \
        --normalisation raw \
        > {log} 2>&1
        """


rule bedgraph_normalised:
    input:
        cooler=rules.merge_counts.output,
    output:
        bedgraph=expand(
            "capcruncher_pileup/04_bedgraph/{{sample}}/norm/{{sample}}.norm.{viewpoint}.bedgraph",
            viewpoint=VIEWPOINT_NAMES,
        ),
    log:
        "logs/bedgraph_norm/{sample}.log",
    params:
        output_prefix=lambda wc, output: f"capcruncher_pileup/04_bedgraph/{wc.sample}/norm/{wc.sample}.norm",
        norm=get_normalisation_from_config(),
        scale_factor=config["normalisation"].get("scale_factor", int(1e6)),
    shell:
        """
        capcruncher \
        interactions \
        pileup \
        {input.cooler} \
        -o {params.output_prefix} \
        {params.norm} \
        --scale-factor {params.scale_factor} \
        > {log} 2>&1
        """


rule bedgraph_to_bigwig:
    input:
        bedgraph="capcruncher_pileup/04_bedgraph/{sample}/{norm}/{sample}.{norm}.{viewpoint}.bedgraph",
    output:
        bigwig="capcruncher_pileup/05_bigwig/{sample}/{norm}/{sample}.{norm}.{viewpoint}.bw",
    log:
        "logs/bedgraph_to_bigwig/{sample}_{norm}_{viewpoint}.log",
    params:
        chrom_sizes=config["genome"]["chrom_sizes"],
    shell:
        """
        sort -k1,1 -k2,2n {input.bedgraph} > {input.bedgraph}.sorted
        bedGraphToBigWig {input.bedgraph}.sorted {params.chrom_sizes} {output.bigwig}
        rm {input.bedgraph}.sorted
        """



