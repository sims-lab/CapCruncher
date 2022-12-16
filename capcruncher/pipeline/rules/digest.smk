rule digest_genome:
    input:
        fasta = config["genome"]["fasta"]
    output:
        bed = "capcruncher_resources/restriction_enzyme_map/genome.digest.bed.gz",
        stats = "capcruncher_statistics/digest_genome/genome_digestion_statistics.txt"
    log:
        "capcruncher_resources/restriction_enzyme_map/genome.digest.log"
    params:
        enzyme_or_site = config["analysis"]["restriction_enzyme"]
    threads:
        4
    shell:
        """
        capcruncher genome digest {input.fasta} -r {params.enzyme_or_site} -o {output.bed}.tmp --sort -l {output.stats} > {log} 2>&1 &&
        pigz -p {threads} {output.bed}.tmp -c > {output.bed} 2> {log}
        """



