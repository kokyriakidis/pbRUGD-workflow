rule svpack_filter_annotated:
    input:
        pbsv_vcf = svpack_input,
        eee_vcf = config['ref']['eee_vcf'],
        gnomadsv_vcf = config['ref']['gnomadsv_vcf'],
        hprc_pbsv_vcf = config['ref']['hprc_pbsv_vcf'],
        gff = config['ref']['ensembl_gff']
    output: f"cohorts/{cohort}/svpack/{cohort}.{ref}.pbsv.svpack.vcf"
    log: f"cohorts/{cohort}/logs/svpack/{cohort}.{ref}.pbsv.svpack.log"
    benchmark: f"cohorts/{cohort}/benchmarks/svpack/{cohort}.{ref}.pbsv.svpack.tsv"
    params:
        min_sv_length = 50
    conda: "envs/svpack.yaml"
    shell:
        """
        (python workflow/scripts/svpack/svpack filter --pass-only \
            --min-svlen {params.min_sv_length} {input.pbsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.eee_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.gnomadsv_vcf} | \
            python workflow/scripts/svpack/svpack match -v - {input.hprc_pbsv_vcf} | \
            python workflow/scripts/svpack/svpack consequence - {input.gff} > {output}) 2> {log}
        """
