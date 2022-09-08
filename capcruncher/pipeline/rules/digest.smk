rule digest_genome:
    input:
        fasta = config["genome"]["fasta"]
    output:
        fasta_digested = "capcruncher_preprocessing/restriction_enzyme_map/genome.digest.bed.gz"
    log:
        "capcruncher_preprocessing/restriction_enzyme_map/genome.digest.log"
    params:
        enzyme_or_site = config["analysis"]["restriction_enzyme"]
    threads:
        4
    shell:
        """
        capcruncher genome digest 
        {input.fasta} 
        -r 
        {params.enzyme_or_site}
        -o
        {output.fasta_digested}.tmp
        --sort
        -l
        {output.stats} 
        &&
        pigz 
        -p 
        {threads} 
        {output.fasta_digested}.tmp
        &&
        mv
        {output.fasta_digested}.tmp
        {output.fasta_digested}
        """



