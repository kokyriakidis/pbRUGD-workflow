rule glnexus:
    input: gvcf_list
    output:
        bcf = temp(f"cohorts/{cohort}/glnexus/{cohort}.{ref}.glnexus.bcf"),
        scratch_dir = temp(f"cohorts/{cohort}/glnexus/{cohort}.{ref}.GLnexus.DB/")
    log: f"cohorts/{cohort}/logs/glnexus/{cohort}.{ref}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/glnexus/{cohort}.{ref}.tsv"
    container: f"docker://quay.io/mlin/glnexus:{config['GLNEXUS_VERSION']}"
    threads: 24
    message: f"Executing {{rule}}: Joint calling variants from {cohort} cohort."
    shell:
        """
        (glnexus_cli --threads {threads} \
            --dir {output.scratch_dir} \
            --config DeepVariant {input} > {output.bcf}) 2> {log}
        """


rule split_glnexus_vcf:
    input:
        vcf = f"cohorts/{cohort}/glnexus/{cohort}.{ref}.glnexus.vcf.gz",
        tbi = f"cohorts/{cohort}/glnexus/{cohort}.{ref}.glnexus.vcf.gz.tbi"
    output: temp(f"cohorts/{cohort}/whatshap/{cohort}.{ref}.regions/{cohort}.{ref}.{{region}}.deepvariant.vcf")
    log: f"cohorts/{cohort}/logs/tabix/query/{cohort}.{ref}.{{region}}.glnexus.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/tabix/query/{cohort}.{ref}.{{region}}.glnexus.vcf.tsv"
    params: region = lambda wildcards: wildcards.region, extra = '-h'
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Extracting {wildcards.region} variants from {input}."
    shell: "tabix {params.extra} {input.vcf} {params.region} > {output} 2> {log}"


rule whatshap_phase:
    input:
        reference = config['ref']['fasta'],
        vcf = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.regions/{cohort}.{ref}.{{chromosome}}.deepvariant.vcf.gz",
        tbi = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.regions/{cohort}.{ref}.{{chromosome}}.deepvariant.vcf.gz.tbi",
        phaseinput = abam_list,
        phaseinputindex = [f"{x}.bai" for x in abam_list]
    output: temp(f"cohorts/{cohort}/whatshap/{cohort}.{ref}.regions/{cohort}.{ref}.{{chromosome}}.deepvariant.phased.vcf.gz")
    log: f"cohorts/{cohort}/logs/whatshap/phase/{cohort}.{ref}.{{chromosome}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/whatshap/phase/{cohort}.{ref}.{{chromosome}}.tsv"
    params:
        chromosome = lambda wildcards: wildcards.chromosome,
        extra = "--indels"
    conda: "envs/whatshap.smk"
    message: "Executing {rule}: Phasing {input.vcf} using {input.phaseinput} for chromosome {wildcards.chromosome}."
    shell:
        """
        (whatshap phase {params.extra} \
            --chromosome {wildcards.chromosome} \
            --output {output} \
            --reference {input.reference} \
            {input.vcf} {input.phaseinput}) > {log} 2>&1
        """

## to include pedigree information in phasing
# extra = """--indels --ped cohorts/{cohort}/{cohort}.ped \
#            --no-genetic-haplotyping --recombination-list cohorts/{cohort}/whatshap/{cohort}.recombination.list"""


rule whatshap_bcftools_concat:
    input:
        calls = expand(f"cohorts/{cohort}/whatshap/{cohort}.{ref}.regions/{cohort}.{ref}.{{region}}.deepvariant.phased.vcf.gz", region=all_chroms),
        indices = expand(f"cohorts/{cohort}/whatshap/{cohort}.{ref}.regions/{cohort}.{ref}.{{region}}.deepvariant.phased.vcf.gz.tbi", region=all_chroms)
    output: f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.phased.vcf.gz"
    log: f"cohorts/{cohort}/logs/bcftools/concat/{cohort}.{ref}.whatshap.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bcftools/concat/{cohort}.{ref}.whatshap.tsv"
    params: "-a -Oz"
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Concatenating WhatsHap phased VCFs: {input.calls}"
    shell: "bcftools concat {params} -o {output} {input.calls} > {log} 2>&1"


rule whatshap_stats:
    input:
        vcf = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.phased.vcf.gz.tbi",
        chr_lengths = config['ref']['chr_lengths']
    output:
        gtf = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.phased.gtf",
        tsv = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.phased.tsv",
        blocklist = f"cohorts/{cohort}//whatshap/{cohort}.{ref}.deepvariant.phased.blocklist"
    log: f"cohorts/{cohort}/logs/whatshap/stats/{cohort}.{ref}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/whatshap/stats/{cohort}.{ref}.tsv"
    conda: "envs/whatshap.yaml"
    message: "Executing {rule}: Calculating phasing stats for {input.vcf}."
    shell:
        """
        (whatshap stats \
            --gtf {output.gtf} \
            --tsv {output.tsv} \
            --block-list {output.blocklist} \
            --chr-lengths {input.chr_lengths} \
            {input.vcf}) > {log} 2>&1
        """


rule link_phased_vcf:
    input:
        vcf = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.phased.vcf.gz.tbi"
    output:
        vcf = f"cohorts/{cohort}/{cohort}.{ref}.deepvariant.phased.vcf.gz",
        tbi = f"cohorts/{cohort}/{cohort}.{ref}.deepvariant.phased.vcf.gz.tbi"
    message: "Linking joint-called, phased vcf to {output.vcf}."
    run:
        to_link = [(input.vcf, output.vcf), (input.tbi, output.tbi)]
        for src, dst in to_link:
            if not os.path.exists(dst):
                os.symlink(("/").join(src.split('/')[2:]), dst)


# TODO
# rule cleanup_whatshap_intermediates:
#     input: f"cohorts/{cohort}/whatshap/{cohort}.{ref}.deepvariant.haplotagged.bam"
#     output: touch(f"cohorts/{cohort}/whatshap/{cohort}.{ref}.removed_intermediates.txt")
#     message: f"Executing {{rule}}: Removing intermediate files for {cohort}."
#     shell:
#         f"""
#         rm -rf cohorts/{cohort}/whatshap/{cohort}.{ref}.regions/
#         """
