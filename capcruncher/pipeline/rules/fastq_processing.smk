import os
import pathlib 

def collate_fastq_files_for_splitting(wc):

    fqs = [fq for fq in sorted(pathlib.Path(".").glob("*.fastq*"))
           if wc.sample in str(fq)]
    return {"fastq_1": fqs[0], "fastq_2": fqs[1]}


checkpoint split_fastqs:
    input:
        **collate_fastq_files(wc)
    output:
        touch("flags/fastq_split/{sample}.sentinel")
    threads:
        2
    params:
        n_reads = str(config["split"].get("n_reads", 1e6)),
        gzip = "--gzip" if COMPRESS_FASTQ else "--no-gzip",
        compression_level = f"--compression_level {config["pipeline"].get('compression', 0)}" if COMPRESS_FASTQ else ""    
    shell:
        """
        capcruncher
        fastq
        split
        {input.fastq_1}
        {input.fastq_2}
        -m
        unix
        -o
        capcruncher_preprocessing/split/{wc.sample}
        -n
        {params.n_reads}
        {params.gzip}
        {params.compression_level}
        -p
        {threads}
        """

def get_split_fastq_names(wc):
    checkpoint_output = checkpoints.split_fastqs.get(**wc).output[0]

    file_names = expand("capcruncher_preprocessing/split/{sample_partition}.fastq.gz", 
                        sample_partition = glob_wildcards(os.path.join(checkpoint_output, "capcruncher_preprocessing/split/{sample_partition}.fastq.gz")))

    return file_names
