def list_svsigs(region):
    """Given a region, return all cohort svsig files from this region."""
    return [f"samples/{sample}/pbsv/svsig/{movie}.{ref}.{region}.svsig.gz" for sample in samples
            for movie in movie_dict[sample]]


localrules: bcftools_concat_pbsv_vcf, cleanup_pbsv_intermediates


rule pbsv_call:
    input:
        svsigs = lambda wildcards: list_svsigs(wildcards.region),
        reference = config['ref']['fasta']
    output: temp(f"cohorts/{cohort}/pbsv/{cohort}.{ref}.chrom_vcfs/{cohort}.{ref}.{{region}}.pbsv.vcf")
    log: f"cohorts/{cohort}/logs/pbsv/call/{cohort}.{ref}.{{region}}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/pbsv/call/{cohort}.{ref}.{{region}}.tsv"
    params:
        region = lambda wildcards: wildcards.region,
        extra = "--ccs -m 20 -A 3 -O 3 -P 20",
        loglevel = "INFO"
    threads: 8
    conda: "envs/pbsv.yaml"
    message: "Executing {rule}: Calling structural variants in {wildcards.region} from SVSIGs: {input.svsigs}"
    shell:
        """
        (pbsv call {params.extra} \
            --log-level {params.loglevel} \
            --num-threads {threads} \
            {input.reference} {input.svsigs} {output}) > {log} 2>&1
        """


rule bcftools_concat_pbsv_vcf:
    input:
        calls = expand(f"cohorts/{cohort}/pbsv/{cohort}.{ref}.chrom_vcfs/{cohort}.{ref}.{{region}}.pbsv.vcf.gz", region=all_chroms),
        indices = expand(f"cohorts/{cohort}/pbsv/{cohort}.{ref}.chrom_vcfs/{cohort}.{ref}.{{region}}.pbsv.vcf.gz.tbi", region=all_chroms)
    output: f"cohorts/{cohort}/pbsv/{cohort}.{ref}.pbsv.vcf"
    log: f"cohorts/{cohort}/logs/bcftools/concat/{cohort}.{ref}.pbsv.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bcftools/concat/{cohort}.{ref}.pbsv.vcf.tsv"
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Concatenating pbsv VCFs: {input.calls}"
    shell: "bcftools concat -a -o {output} {input.calls}"


rule link_pbsv_vcf:
    input:
        vcf = f"cohorts/{cohort}/pbsv/{cohort}.{ref}.pbsv.vcf.gz",
        tbi = f"cohorts/{cohort}/pbsv/{cohort}.{ref}.pbsv.vcf.gz.tbi"
    output:
        vcf = f"cohorts/{cohort}/{cohort}.{ref}.pbsv.vcf.gz",
        tbi = f"cohorts/{cohort}/{cohort}.{ref}.pbsv.vcf.gz.tbi"
    message: "Linking joint-called vcf to {output.vcf}."
    run:
        for src, dst in to_link:
            if not os.path.exists(dst):
                os.symlink(("/").join(src.split('/')[3:]), dst)


# TODO
# rule cleanup_pbsv_intermediates:
#     input: f"cohorts/{cohort}/{cohort}.{ref}.pbsv.vcf.gz"
#     output: touch(f"cohorts/{cohort}/pbsv/{cohort}.{ref}.removed_intermediates.txt")
#     message: f"Executing {{rule}}: Removing intermediate folder for {cohort}."
#     shell: f"rm -rf cohorts/{cohort}/pbsv/{cohort}.{ref}.chrom_vcfs/"
