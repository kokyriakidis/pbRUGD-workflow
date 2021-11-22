localrules: reformat_ensembl_gff, generate_lof_lookup, generate_clinvar_lookup


rule reformat_ensembl_gff:
    output: config['ref']['ensembl_gff']
    log: "logs/reformat_ensemble_gff.log"
    params: url = config['ref']['ensembl_gff_url']
    conda: "envs/htslib.yaml"
    message: "Executing {rule}: Downloaded and reformatting ensembl GFF to {output}."
    shell:
        """
        (wget -qO - {params.url} | zcat \
            | awk -v OFS="\t" '{{ if ($1=="##sequence-region") && ($2~/^G|K/) {{ print $0; }} else if ($0!~/G|K/) {{ print "chr" $0; }} }}' \
            | bgzip > {output}) > {log} 2>&1
        """


rule generate_lof_lookup:
    output: config['lof_lookup']
    log: "logs/generate_lof_lookup.log"
    params: url = config['lof_lookup_url']
    message: "Executing {rule}: Generating a lookup table for loss-of-function metrics at {output}."
    shell:
        """
        (wget -qO - {params.url} | zcat | cut -f 1,21,24 | tail -n+2 \
            | awk "{{ printf(\\"%s\\tpLI=%.3g;oe_lof=%.5g\\n\\", \$1, \$2, \$3) }}" > {output}) > {log} 2>&1
        """


rule generate_clinvar_lookup:
    output: config['clinvar_lookup']
    log: "logs/generate_clinvar_lookup.log"
    params: url = config['clinvar_lookup_url']
    message: "Executing {rule}: Generating a lookup table for clinvar gene descriptions at {output}."
    shell: "(wget -qO - {params.url} | cut -f 2,5 | grep -v ^$'\t' > {output}) > {log} 2>&1"


rule bcftools_norm:
    input:
        vcf = slivar_input,
        tbi = slivar_input + ".tbi"
    output: temp(f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.norm.bcf")
    log: f"cohorts/{cohort}/logs/bcftools/norm/{cohort}.{ref}.deepvariant.phased.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/bcftools/norm/{cohort}.{ref}.deepvariant.phased.vcf.tsv"
    params: f"--multiallelics - --output-type b --fasta-ref {config['ref']['fasta']}"
    conda: "envs/bcftools.yaml"
    message: "Executing {rule}: Splitting multiallelic sites and normalizing indels for {input.vcf}."
    shell: "(bcftools norm {params} {input.vcf} -o {output}) > {log} 2>&1"


if singleton:
    # singleton
    slivar_filters = [
        "--info 'variant.FILTER==\"PASS\" && INFO.gnomad_af < 0.01 && INFO.hprc_af < 0.01 && INFO.gnomad_nhomalt < 5 && INFO.hprc_nhomalt < 5'",
        "--family-expr 'recessive:fam.every(segregating_recessive)'",
        "--family-expr 'x_recessive:(variant.CHROM == \"chrX\") && fam.every(segregating_recessive_x)'",
        "--family-expr 'dominant:fam.every(segregating_dominant) && INFO.gnomad_ac < 5 && INFO.hprc_ac < 5'",
        "--family-expr 'x_dominant:(variant.CHROM == \"chrX\") && fam.every(segregating_dominant_x) && INFO.gnomad_ac < 5 && INFO.hprc_ac < 5'",
        "--sample-expr 'comphet_side:sample.het && sample.GQ > 5'"
        ]
else:
    # trio cohort
    slivar_filters = [
        "--info 'variant.FILTER==\"PASS\" && INFO.gnomad_af < 0.01 && INFO.hprc_af < 0.01 && INFO.gnomad_nhomalt < 5 && INFO.hprc_nhomalt < 5'",
        "--family-expr 'recessive:fam.every(segregating_recessive)'",
        "--family-expr 'x_recessive:(variant.CHROM == \"chrX\") && fam.every(segregating_recessive_x)'",
        "--family-expr 'dominant:fam.every(segregating_dominant) && INFO.gnomad_ac < 5 && INFO.hprc_ac < 5'",
        "--family-expr 'x_dominant:(variant.CHROM == \"chrX\") && fam.every(segregating_dominant_x) && INFO.gnomad_ac < 5 && INFO.hprc_ac < 5'",
        "--trio 'comphet_side:comphet_side(kid, mom, dad) && kid.affected'"
    ]


rule slivar_small_variant:
    input:
        bcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.norm.bcf",
        csi = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.norm.bcf.csi",
        ped = f"cohorts/{cohort}/{cohort}.ped",
        gnomad_af = {config['ref']['gnomad_gnotate']},
        hprc_af = {config['ref']['hprc_dv_gnotate']},
        js = config['slivar_js'],
        gff = config['ref']['ensembl_gff'],
        ref = config['ref']['fasta']
    output: f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf"
    log: f"cohorts/{cohort}/logs/slivar/filter/{cohort}.{ref}.deepvariant.phased.slivar.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/filter/{cohort}.{ref}.deepvariant.phased.slivar.tsv"
    params: filters = slivar_filters
    threads: 12
    conda: "envs/slivar.yaml"
    message: "Executing {rule}: Annotating {input.bcf} and applying filters."
    shell:
        """
        (pslivar --processes {threads} \
            --fasta {input.ref}\
            --pass-only \
            --js {input.js} \
            {params.filters} \
            --gnotate {input.gnomad_af} \
            --gnotate {input.hprc_af} \
            --vcf {input.bcf} \
            --ped {input.ped} \
            | bcftools csq -l -s - --ncsq 40 \
                -g {input.gff} -f {input.ref} - -o {output}) > {log} 2>&1
        """


skip_list = [
    'non_coding_transcript',
    'intron',
    'non_coding_transcript',
    'non_coding',
    'upstream_gene',
    'downstream_gene',
    'non_coding_transcript_exon',
    'NMD_transcript',
    '5_prime_UTR',
    '3_prime_UTR'
    ]


rule slivar_compound_hets:
    input: 
        vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz",
        tbi = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz.tbi",
        ped = f"cohorts/{cohort}/{cohort}.ped"
    output:
        vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf"
    log: f"cohorts/{cohort}/logs/slivar/compound-hets/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/compound-hets/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf.tsv"
    params: skip = ",".join(skip_list)
    conda: "envs/slivar.yaml"
    message: f"Executing {{rule}}: Finding compound hets in {cohort}."
    shell:
        """
        (slivar compound-hets \
            --skip {params.skip} \
            --vcf {input.vcf} \
            --sample-field comphet_side \
            --ped {input.ped} \
            --allow-non-trios \
            | python3 workflow/scripts/add_comphet_phase.py \
            > {output.vcf}) > {log} 2>&1
        """


info_fields = [
    'gnomad_af',
    'hprc_af',
    'gnomad_nhomalt',
    'hprc_nhomalt',
    'gnomad_ac',
    'hprc_ac'
    ]


rule slivar_tsv:
    input:
        filt_vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.vcf.gz",
        comphet_vcf = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.vcf.gz",
        ped = f"cohorts/{cohort}/{cohort}.ped",
        lof_lookup = config['lof_lookup'],
        clinvar_lookup = config['clinvar_lookup'],
        phrank_lookup = f"cohorts/{cohort}/{cohort}_phrank.tsv"
    output:
        filt_tsv = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.tsv",
        comphet_tsv = f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.phased.slivar.compound-hets.tsv"
    log: f"cohorts/{cohort}/logs/slivar/tsv/{cohort}.{ref}.log"
    benchmark: f"cohorts/{cohort}/benchmarks/slivar/tsv/{cohort}.{ref}.tsv"
    params: info = "".join([f"--info-field {x} " for x in info_fields])
    conda: "envs/slivar.yaml"
    message: "Executing {rule}: Converting annotated VCFs to TSVs for easier interpretation."
    shell:
        """
        (slivar tsv \
            {params.info} \
            --sample-field dominant \
            --sample-field x_dominant \
            --sample-field recessive \
            --sample-field x_recessive \
            --csq-field BCSQ \
            --gene-description {input.lof_lookup} \
            --gene-description {input.clinvar_lookup} \
            --gene-description {input.phrank_lookup} \
            --ped {input.ped} \
            --out /dev/stdout \
            {input.filt_vcf} \
            | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
            > {output.filt_tsv}
        slivar tsv \
            {params.info} \
            --sample-field slivar_comphet \
            --info-field slivar_comphet \
            --csq-field BCSQ \
            --gene-description {input.lof_lookup} \
            --gene-description {input.clinvar_lookup} \
            --gene-description {input.phrank_lookup} \
            --ped {input.ped} \
            --out /dev/stdout \
            {input.comphet_vcf} \
            | sed '1 s/gene_description_1/lof/;s/gene_description_2/clinvar/;s/gene_description_3/phrank/;' \
            > {output.comphet_tsv}) > {log} 2>&1
        """
