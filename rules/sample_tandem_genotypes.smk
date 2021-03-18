rule last_align:
    input:
        bam = f"samples/{sample}/whatshap/{sample}.{ref}.deepvariant.haplotagged.bam",
        bed = config['tg_targets'],
        score_matrix = config['score_matrix']
    output: temp(f"samples/{sample}/tandem-genotypes/{sample}.maf.gz")
    log: f"samples/{sample}/logs/last/align/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/last/align/{sample}.tsv"
    conda: "envs/last.yaml"
    params: 
        last_index = config['ref']['last_index'],
        extra = "-C2"
    threads: 24
    message: "Executing {rule}: Aligning {input.query} to {input.last_index} using lastal with {input.score_matrix} score matrix."
    shell:
        """
        (samtools view -@7 -bL {input.bed} {input.bam} | samtools fasta \
         | lastal -P16 -p {input.score_matrix} {params.extra} \
               {params.last_index} {input.query} \
               | last-split | bgzip > {output}) > {log} 2>&1
        """


rule tandem_genotypes:
    input:
        maf: f"samples/{sample}/tandem-genotypes/{sample}.maf.gz",
        repeats: config['tandem_repeats_targets']
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem_repeats.txt"
    log: f"samples/{sample}/logs/tandem_genotypes/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/tandem_genotypes/{sample}.tsv"
    conda: "envs/tandem_repeats.yaml"
    message: "Executing {rule}: Genotyping tandem repeats from {input.repeats} regions in {input.maf}."
    shell: "tandem-genotypes {input.repeats} {input.maf} > {output}"


rule tandem_genotypes_plot:
    input: f"samples/{sample}/tandem-genotypes/{sample}.tandem_repeats.txt"
    output: f"samples/{sample}/tandem-genotypes/{sample}.tandem_repeats.pdf"
    log: f"samples/{sample}/logs/tandem_genotypes/plot/{sample}.log"
    benchmark: f"samples/{sample}/benchmarks/tandem_genotypes/plot/{sample}.tsv"
    conda: "envs/tandem_repeats.yaml"
    params: top_N_plots = 50
    message: "Executing {rule}: Plotting tandem repeats from {input}."
    shell: "tandem-genotypes-plot -n {params.top_N_plots} {input} {output}"
