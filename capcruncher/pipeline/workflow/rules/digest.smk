rule digest_genome:
    input:
        fasta=config["genome"]["fasta"],
    output:
        bed="capcruncher_output/resources/restriction_fragments/genome.digest.bed.gz",
    log:
        "capcruncher_output/resources/restriction_fragments/genome.digest.log",
    params:
        enzyme_or_site=config["analysis"]["restriction_enzyme"],
    threads: 4
    resources:
        mem_mb=2000,
    shell:
        """
        capcruncher genome digest {input.fasta} -r {params.enzyme_or_site} -o {output.bed}.tmp --sort > {log} 2>&1 &&
        pigz -p {threads} {output.bed}.tmp -c > {output.bed} 2> {log}
        rm {output.bed}.tmp
        """
