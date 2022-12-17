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


rule done:
    input:
        "capcruncher_pileup/03_counts/{sample}.hdf5",
    output:
        touch("{sample}.done"),


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
        reporters \
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
        reporters \
        store \
        bins \
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
    threads: 4
    shell:
        """
        capcruncher \
        reporters \
        merge \
        {input} \
        -o {output} \
        -p {threads} \
        > {log} 2>&1
        """
