import pathlib


rule fastq_rename:
    input:
        fq1=lambda wc: FASTQ_SAMPLES.translation[f"{wc.sample}_1.fastq.gz"],
        fq2=lambda wc: FASTQ_SAMPLES.translation[f"{wc.sample}_2.fastq.gz"],
    output:
        fq1="capcruncher_output/interim/fastq/{sample}_1.fastq.gz",
        fq2="capcruncher_output/interim/fastq/{sample}_2.fastq.gz",
    log:
        "capcruncher_output/logs/fastq_rename/{sample}.log",
    shell:
        """
        ln -s $(realpath {input.fq1}) {output.fq1} \
            && ln -s $(realpath {input.fq2}) {output.fq2}
        """


checkpoint split:
    input:
        fq1=rules.fastq_rename.output.fq1,
        fq2=rules.fastq_rename.output.fq2,
    output:
        directory("capcruncher_output/interim/fastq/split/{sample}"),
    log:
        "capcruncher_output/logs/split/{sample}.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(1, attempt),
        runtime=lambda wildcards, attempt: scale_resource(180, attempt),
    params:
        prefix="capcruncher_output/interim/fastq/split/{sample}/{sample}",
        n_reads=str(config["split"].get("n_reads", 1e6)),
        method=default_fastq_split_method(),
    shell:
        """
        mkdir {output} \
            && capcruncher \
                fastq \
                split \
                {input.fq1} \
                {input.fq2} \
                -m \
                {params.method} \
                -o \
                {params.prefix} \
                -n \
                {params.n_reads} \
                --gzip \
                -p \
                {threads} \
                >{log} 2>&1
        """


checkpoint deduplication:
    input:
        fq1=lambda wc: get_fastq_split_1(wc)["fq1"],
        fq2=lambda wc: get_fastq_split_1(wc)["fq2"],
    output:
        fastq_dir=directory("capcruncher_output/interim/fastq/deduplicated/{sample}/"),
        stats="capcruncher_output/interim/statistics/deduplication/data/{sample}.deduplication.json",
    log:
        "capcruncher_output/logs/deduplication_fastq/{sample}.log",
    threads: max(1, workflow.cores // 2)
    resources:
        mem=lambda wildcards, attempt: scale_memory(2, attempt),
    params:
        prefix_fastq="capcruncher_output/interim/fastq/deduplicated/{sample}/",
        fq1_args=lambda wc, input: " ".join(f"-1 {f}" for f in (input.fq1 if isinstance(input.fq1, list) else [input.fq1])),
        fq2_args=lambda wc, input: " ".join(f"-2 {f}" for f in (input.fq2 if isinstance(input.fq2, list) else [input.fq2])),
    shell:
        """
        mkdir -p {params.prefix_fastq} \
            && capcruncher fastq deduplicate {params.fq1_args} {params.fq2_args} -o {params.prefix_fastq} --statistics {output.stats} --sample-name {wildcards.sample} >{log} 2>&1
        """


rule trim:
    input:
        fq1=lambda wc: get_deduplicated_fastq_pair(wc)["fq1"],
        fq2=lambda wc: get_deduplicated_fastq_pair(wc)["fq2"],
    output:
        trimmed1=temp(
            "capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_1.fastq.gz"
        ),
        trimmed2=temp(
            "capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_2.fastq.gz"
        ),
    log:
        "capcruncher_output/logs/trimming/{sample}_{part}.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(2, attempt),
    params:
        outdir="capcruncher_output/interim/fastq/trimmed/{sample}/",
    shell:
        """
        trim_galore --cores {threads} --trim-n --paired --output_dir {params.outdir} {input.fq1} {input.fq2} >>{log} 2>&1 \
            && mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_1_val_1.fq.gz {output.trimmed1} \
            && mv {params.outdir}/{wildcards.sample}_part{wildcards.part}_2_val_2.fq.gz {output.trimmed2}
        """


rule flash:
    input:
        fq1="capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_1.fastq.gz",
        fq2="capcruncher_output/interim/fastq/trimmed/{sample}/{sample}_part{part}_2.fastq.gz",
    output:
        flashed="capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.extendedFrags.fastq.gz",
        pe1="capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_1.fastq.gz",
        pe2="capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.notCombined_2.fastq.gz",
        hist=temp(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.hist"
        ),
        histogram=temp(
            "capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}.histogram"
        ),
    log:
        "capcruncher_output/logs/flash/{sample}_{part}.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(1, attempt),
    params:
        outdir="capcruncher_output/interim/fastq/flashed/{sample}/{sample}_part{part}",
    shell:
        """
        flash2 {input.fq1} {input.fq2} -o {params.outdir} -t {threads} -z >{log} 2>&1
        """


checkpoint rebalance_partitions_combined:
    input:
        flashed=lambda wc: get_flashed_fastq(wc),
    output:
        fastq_dir=directory(
            "capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/"
        ),
        sentinel=touch(
            "capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/.complete.sentinel"
        ),
    log:
        "capcruncher_output/logs/rebalance_partitions/{sample}_flashed.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(1, attempt),
    params:
        prefix=lambda wildcards, output: pathlib.Path(output.fastq_dir)
        / wildcards.sample,
        suffix=lambda wc: f"_flashed",
        fq=lambda wc: ",".join(get_flashed_fastq(wc)),
        n_reads=str(config["split"].get("n_reads", 1e6)),
    shell:
        """
        mkdir -p {output.fastq_dir} \
            && capcruncher \
                fastq \
                split \
                {params.fq} \
                -m \
                unix \
                -o \
                {params.prefix} \
                -n \
                {params.n_reads} \
                --gzip \
                -p \
                {threads} \
                --suffix \
                {params.suffix} \
                >{log} 2>&1 \
            && touch {output.sentinel}
        """


checkpoint rebalance_partitions_pe:
    input:
        fq=get_pe_fastq,
    output:
        fastq_dir=directory("capcruncher_output/interim/fastq/rebalanced/{sample}/pe"),
        sentinel=touch(
            "capcruncher_output/interim/fastq/rebalanced/{sample}/pe/.complete.sentinel"
        ),
    log:
        "capcruncher_output/logs/rebalance_partitions/{sample}_pe.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(1, attempt),
    params:
        prefix=lambda wildcards, output: pathlib.Path(output.fastq_dir)
        / wildcards.sample,
        suffix=lambda wc: f"_pe",
        n_reads=str((config["split"].get("n_reads", 1e6) // 2)),
        fq1=lambda wc: ",".join(separate_pe_fastq(wc)[1]),
        fq2=lambda wc: ",".join(separate_pe_fastq(wc)[2]),
    shell:
        """
        mkdir -p {output.fastq_dir} \
            && capcruncher \
                fastq \
                split \
                {params.fq1} \
                {params.fq2} \
                -m \
                unix \
                -o \
                {params.prefix} \
                -n \
                {params.n_reads} \
                --gzip \
                -p \
                {threads} \
                --suffix \
                {params.suffix} \
                >{log} 2>&1 \
            && touch {output.sentinel}
        """


rule digest_flashed_combined:
    input:
        flashed="capcruncher_output/interim/fastq/rebalanced/{sample}/flashed/{sample}_part{part}_flashed_1.fastq.gz",
    output:
        digested=temp(
            "capcruncher_output/interim/fastq/digested/{sample}/{sample}_part{part}_flashed.fastq.gz"
        ),
        statistics="capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_flashed.json",
    log:
        "capcruncher_output/logs/digestion/{sample}_{part}.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(2, attempt),
    params:
        restriction_site=config["analysis"]["restriction_enzyme"],
    shell:
        """
        capcruncher \
            fastq \
            digest \
            {input.flashed} \
            -o \
            {output.digested} \
            -m \
            flashed \
            -r \
            {params.restriction_site} \
            --minimum_slice_length \
            18 \
            --statistics \
            {output.statistics} \
            --sample-name \
            {wildcards.sample} \
            >{log} 2>&1
        """


rule digest_flashed_pe:
    input:
        pe1="capcruncher_output/interim/fastq/rebalanced/{sample}/pe/{sample}_part{part}_pe_1.fastq.gz",
        pe2="capcruncher_output/interim/fastq/rebalanced/{sample}/pe/{sample}_part{part}_pe_2.fastq.gz",
    output:
        digested=temp(
            "capcruncher_output/interim/fastq/digested/{sample}/{sample}_part{part}_pe.fastq.gz"
        ),
        statistics="capcruncher_output/interim/statistics/digestion/data/{sample}_part{part}_pe.json",
    log:
        "capcruncher_output/logs/digestion/{sample}_{part}.log",
    threads: 4
    resources:
        mem=lambda wildcards, attempt: scale_memory(2, attempt),
    params:
        restriction_site=config["analysis"]["restriction_enzyme"],
    shell:
        """
        capcruncher \
            fastq \
            digest \
            {input.pe1} \
            {input.pe2} \
            -o \
            {output.digested} \
            -m \
            pe \
            -r \
            {params.restriction_site} \
            --minimum_slice_length \
            18 \
            --statistics \
            {output.statistics} \
            --sample-name \
            {wildcards.sample} \
            >{log} 2>&1
        """


localrules:
    fastq_rename,
