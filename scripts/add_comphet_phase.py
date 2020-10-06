#!/usr/bin/env python3
"""
For slivar_comphet tags, note whether pairs are cis/trans/unknown
"""

__author__ = "William Rowell"
__version__ = "0.1.0"
__license__ = "MIT"


import vcfpy
from collections import defaultdict


def import_records(reader):
    """Import VCF and return records list and lookup.

    Given a vcfpy.Reader object, return a list of records in the original order
    as well as a dictionary with unique variant identifiers as key and phase
    set and genotype as values.
    """
    records = list(reader)
    lookup = defaultdict(dict)
    for record in records:
        varkey = (record.CHROM, str(record.POS),
                  record.REF, record.ALT[0].value)
        for sample in record.calls:
            samplename = sample.sample
            GT = sample.data['GT']
            PS = '0'
            if sample.is_phased:
                PS = str(sample.data['PS'])
        lookup[varkey][samplename] = (PS, GT)
    return records, lookup


def compare_phase(slivar_comphet, calls, lookup, phase_table):
    """Return whether slivar_compet and this variant are on same phase.

    Given current slivar_comphet record and variant calls for current record
    return cis if variants are on same haplotype, trans if on opposite
    haplotypes, and unknown otherwise."""
    sample, gene, chid, chrom, pos, ref, alt = slivar_comphet.split('/')
    ch_PS, ch_GT = lookup[(chrom, pos, ref, alt)][sample]
    # look up the phase set and genotype of this variant
    if not calls[sample].is_phased:
        phase_table["_".join(["slivar_comphet", str(chid)])] = 'unknown'
        return 'unknown'
    this_GT = calls[sample].data['GT']
    this_PS = str(calls[sample].data['PS'])
    if ch_PS == this_PS:
        if ch_GT == this_GT:
            phase_table["_".join(["slivar_comphet", str(chid)])] = 'cis'
            return 'cis'
        elif ch_GT != this_GT:
            phase_table["_".join(["slivar_comphet", str(chid)])] = 'trans'
            return 'trans'
    else:
        phase_table["_".join(["slivar_comphet", str(chid)])] = 'unknown'
        return 'unknown'


def main():
    # read from stdin
    reader = vcfpy.Reader.from_path('/dev/stdin')
    records, lookup = import_records(reader)

    phase_table = {}
    # write to stdout
    with vcfpy.Writer.from_path('/dev/stdout', reader.header) as writer:
        for record in records:
            calls = {x.sample: x for x in record.calls}
            for ix, slivar_comphet in enumerate(record.INFO['slivar_comphet']):
                phase = compare_phase(slivar_comphet, calls,
                                      lookup, phase_table)
                record.INFO['slivar_comphet'][ix] = \
                    "/".join([slivar_comphet, phase])
            writer.write_record(record)


if __name__ == "__main__":
    """ This is executed when run from the command line """
    main()
