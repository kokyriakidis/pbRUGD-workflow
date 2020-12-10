rule deepvariant_round1:
    input:
        bams = abams,
        bais = [f"{x}.bai" for x in abams],
        reference = config['ref']['fasta']
    output:
        vcf = f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.vcf.gz",
        vcf_index = f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.vcf.gz.tbi",
        report = f"samples/{sample}/deepvariant_intermediate/{sample}.{ref}.deepvariant.visual_report.html"
    log: f"samples/{sample}/logs/deepvariant_intermediate/{sample}.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/deepvariant_intermediate/{sample}.{ref}.tsv"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params: reads = ','.join(abams)
    threads: 60
    message: "Executing {rule}: DeepVariant round1 for {input.bams}."
    shell:
        """
        (/opt/deepvariant/bin/run_deepvariant \
            --model_type PACBIO \
            --num_shards {threads} \
            --ref {input.reference} \
            --reads {params.reads} \
            --output_vcf {output.vcf}) > {log} 2>&1
        """


haplotagged_abams = [f"samples/{sample}/whatshap_intermediate/{sample}.{ref}.{movie}.deepvariant.haplotagged.bam" for movie in movies]


rule deepvariant_round2:
    input:
        bams = haplotagged_abams,
        bais = [f"{x}.bai" for x in haplotagged_abams],
        reference = config['ref']['fasta']
    output:
        vcf = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz",
        vcf_index = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz.tbi",
        gvcf = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz",
        gvcf_index = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz.tbi",
        report = f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.visual_report.html"
    log: f"samples/{sample}/logs/deepvariant/{sample}.{ref}.log"
    benchmark: f"samples/{sample}/benchmarks/deepvariant/{sample}.{ref}.tsv"
    container: f"docker://google/deepvariant:{config['DEEPVARIANT_VERSION']}"
    params: reads = ','.join(haplotagged_abams)
    threads: 60
    message: "Executing {rule}: DeepVariant round1 for {input.bams}."
    shell:
        """
        (/opt/deepvariant/bin/run_deepvariant \
            --model_type PACBIO \
            --make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true" \
            --num_shards {threads} \
            --ref {input.reference} \
            --reads {params.reads} \
            --output_vcf {output.vcf} \
            --output_gvcf {output.gvcf}) > {log} 2>&1
        """
# command above is for DV v1.0
# if using DV v1.1,
#   replace `--make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true" \`
#   with `--use_hp_information \`


rule bcftools_stats:
    input: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.gz"
    output: f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.vcf.stats.txt"
    log: f"samples/{sample}/logs/bcftools/stats/{sample}.{ref}.deepvariant.vcf.log"
    benchmark: f"samples/{sample}/benchmarks/bcftools/stats/{sample}.{ref}.deepvariant.vcf.tsv"
    params: f"--fasta-ref {config['ref']['fasta']} --apply-filters PASS -s {sample}"
    threads: 4
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Calculating VCF statistics for {input}."
    shell: "(bcftools stats --threads 3 {params} {input} > {output}) > {log} 2>&1"
