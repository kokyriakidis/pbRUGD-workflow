import os
import re
from pathlib import Path
from collections import defaultdict


shell.prefix("set -o pipefail; umask 002; ")           # set g+w
configfile: "workflow/reference.yaml"                  # reference information
configfile: "workflow/config.yaml"                     # general configuration
configfile: "100humans-data/cohorts/all_as_dict.yaml"  # list of all cohorts and samples TODO

all_chroms = config['ref']['autosomes'] + config['ref']['sex_chrom'] + config['ref']['mit_chrom']

# cohort will be provided at command line with `--config cohort=$COHORT`
cohort = config['cohort']
ref = config['ref']['shortname']
samples = []
print(f"Processing cohort {cohort} with reference {ref}.")

# find all samples in cohort
for status in ("affecteds", "unaffecteds"):
    if status in config[cohort]:
        samples.extend([s['id'] for s in config[cohort][status]])
singleton = False
if len(samples) == 0:
    print(f"No samples in {cohort}.")
elif len(samples) == 1:
    singleton = True
    # for singletons, use the VCFs from the sample folder
    # link sample pbsv and deepvariant vcfs into cohort folder
    to_link = [
        (f"samples/{samples[0]}/pbsv/{samples[0]}.{ref}.pbsv.vcf.gz",
         f"cohorts/{cohort}/{cohort}.{ref}.pbsv.vcf.gz"),
        (f"samples/{samples[0]}/pbsv/{samples[0]}.{ref}.pbsv.vcf.gz.tbi",
         f"cohorts/{cohort}/{cohort}.{ref}.pbsv.vcf.gz.tbi"),
        (f"samples/{samples[0]}/whatshap/{samples[0]}.{ref}.deepvariant.phased.vcf.gz",
         f"cohorts/{cohort}/{cohort}.{ref}.deepvariant.phased.vcf.gz"),
        (f"samples/{samples[0]}/whatshap/{samples[0]}.{ref}.deepvariant.phased.vcf.gz.tbi",
         f"cohorts/{cohort}/{cohort}.{ref}.deepvariant.phased.vcf.gz.tbi")
    ]
    for src, dst in to_link:
        if not os.path.exists(dst):
            os.symlink("../../" + src, dst)
print(f"Samples in cohort: {samples}.")

# scan samples/*/aligned to generate a dict-of-lists-of-movies for 
pattern = re.compile(r'samples/(?P<sample>[A-Za-z0-9_-]+)/aligned/(?P<movie>m\d{5}[U]?_\d{6}_\d{6})\.(?P<reference>.*).bam')
movie_dict = defaultdict(list)
abam_list = []
for infile in Path(f"samples").glob('**/aligned/*.bam'):
    match = pattern.search(str(infile))
    if match and (match.group('sample') in samples) and (match.group('reference') == ref):
        movie_dict[match.group('sample')].append([match.group('movie')])
        abam_list.append(infile)
gvcf_list = [f"samples/{sample}/deepvariant/{sample}.{ref}.deepvariant.g.vcf.gz" for sample in samples]

# build a list of targets
targets = []
include: 'rules/cohort_common.smk'

# generate a cohort level pbsv vcf
include: 'rules/cohort_pbsv.smk'
targets.extend([f"cohorts/{cohort}/{cohort}.{ref}.pbsv.{suffix}"
                for suffix in ['vcf.gz', 'vcf.gz.tbi']])

# TODO: annotate and filter pbsv vcf

# generate a cohort level deepvariant vcf
include: 'rules/cohort_glnexus.smk'
targets.extend([f"cohorts/{cohort}/{cohort}.{ref}.deepvariant.phased.{suffix}"
                for suffix in ['vcf.gz', 'vcf.gz.tbi']])

# annotate and filter deepvariant vcf
include: 'rules/cohort_slivar.smk'
targets.extend([f"cohorts/{cohort}/slivar/{cohort}.{ref}.deepvariant.{infix}.{suffix}"
                for infix in ['slivar', 'slivar.compound-hets']
                for suffix in ['vcf.gz', 'vcf.gz.tbi', 'tsv']])


localrules: all, md5sum


rule all:
    input: targets + [f"{x}.md5" for x in targets]


rule md5sum:
    input: "{prefix}"
    output: "{prefix}.md5"
    message: "Creating md5 checksum for {input}."
    shell: "md5sum {input} > {output}"


