# pbRUGD-workflow

## Workflow for the comprehensive detection and prioritization of variants in human genomes with PacBio HiFi reads

_Note_: Code coming to this repo soon.  Currently removing last remaining internal dependencies.

## Authors

- William Rowell ([@williamrowell](https://github.com/williamrowell))
- Aaron Wenger ([@amwenger](https://github.com/amwenger))

## Description

This repo consists of three [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflows:

1. [process_smrtcells](#process_smrtcells)
2. [process_samples](#process_samples)
3. [process_cohorts](#process_cohorts)

### `process_smrtcells`

- find new HiFi BAMs or FASTQs under `smrtcells/ready/*/`
- align HiFi reads to reference (GRCh38 by default) with [pbmm2](https://github.com/PacificBiosciences/pbmm2)
- calculate aligned coverage depth with [mosdepth](https://github.com/brentp/mosdepth)
- calculate depth ratios (chrX:chrY, chrX:chr2) from mosdepth summary to check for sample swaps
- calculate depth ratio (chrM:chr2) from mosdepth summary to check for consistency between runs
- count kmers in HiFi reads using [jellyfish](https://github.com/gmarcais/Jellyfish), dump and export modimers

### `process_sample`

- launch once sample has been sequenced to sufficient depth
- discover and call structural variants with [pbsv](https://github.com/PacificBiosciences/pbsv)
- call small variants with [DeepVariant](https://github.com/google/deepvariant)
- phase small variants with [WhatsHap](https://github.com/whatshap/whatshap/)
- merge per SMRT Cell BAMs and tag merged bam with haplotype based on WhatsHap phased DeepVariant variant calls
- merge jellyfish kmer counts
- assemble reads with [hifiasm](https://github.com/chhylp123/hifiasm)
- align assembly to reference with [minimap2](https://github.com/lh3/minimap2)
- check for sample swaps by calculate consistency of kmers between sequencing runs

### `process_cohort`

- launched once all samples in cohort have been processed
- if singleton
  - link sample level structural variant VCF into cohort folder
  - link phased small variant VCF into cohort folder
- if cohort
  - jointly call structural variants with pbsv
  - jointly call small variants with [GLnexus](https://github.com/dnanexus-rnd/GLnexus)
- using [slivar](https://github.com/brentp/slivar)
  - annotate variant calls with population frequency from [gnomAD](https://gnomad.broadinstitute.org) and HPRC variant databases
  - filter variant calls according to population frequency and inheritance patterns
  - detect possible compound heterozygotes, and filter to remove cis-combinations

## Dependencies

- conda
- singularity >= 3.5.3
