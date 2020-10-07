#!/usr/bin/env python3
"""
Create PED file for cohortid from 100Humans cohortyaml
"""

__author__ = "William Rowell"
__version__ = "0.1.0"
__license__ = "MIT"

import argparse
import yaml
import csv


"""
Family_ID = cohortid
Individual_ID

# assumed that maternal_id is listed first in parents
Paternal_ID - '.' for unknown
Maternal_ID - '.' for unknown

Sex - '1'=male; '2'=female; '.''=unknown
Phenotype - '1'=unaffected, '2'=affected
"""

FIELDNAMES = ['Family_ID', 'Individual_ID',
              'Paternal_ID', 'Maternal_ID',
              'Sex', 'Phenotype']
SEXES = {'MALE': 1, 'FEMALE': 2}
PHENOTYPES = {'unaffecteds': 1, 'affecteds': 2}


def find_cohort(args):
    """Find entry for cohortid within cohortyaml."""
    with open(args.cohortyaml, 'r') as yamlfile:
        cohort_list = yaml.load(yamlfile, Loader=yaml.FullLoader)
    for cohort in cohort_list:
        if cohort['id'] == args.cohortid:
            return cohort
    print(f"Cohort {args.cohortid} not found in {args.cohortyaml}.") and exit


def parse_yaml(args):
    """Parse cohortyaml, returning list of dicts for each individual in cohortid."""
    rows = []
    cohort = find_cohort(args)
    for phenotype in PHENOTYPES.keys():
        if phenotype in cohort:
            for individual in range(len(cohort[phenotype])):
                # if sex is absent, set to unknown value
                if 'sex' in cohort[phenotype][individual]:
                    ind_sex = SEXES[cohort[phenotype][individual]['sex']]
                else:
                    ind_sex = '.'
                # assume that maternal_id is always listed first in parents
                if ('parents' in cohort[phenotype][individual]) \
                        and (len(cohort[phenotype][individual]['parents']) == 2):
                    mat_id = cohort[phenotype][individual]['parents'][0]
                    pat_id = cohort[phenotype][individual]['parents'][1]
                else:
                    mat_id = '.'
                    pat_id = '.'
                row = {
                    'Family_ID': args.cohortid,
                    'Individual_ID': cohort[phenotype][individual]['id'],
                    'Paternal_ID': pat_id,
                    'Maternal_ID': mat_id,
                    'Sex': ind_sex,
                    'Phenotype': PHENOTYPES[phenotype]
                }
                rows.append(row)
    return rows


def write_ped(args, rows):
    """Write list of dicts from parse_yaml."""
    with open(args.pedigree, 'w') as pedigreefile:
        writer = csv.DictWriter(pedigreefile, fieldnames=FIELDNAMES,
                                delimiter='\t')
        for row in rows:
            writer.writerow(row)


def main(args):
    """"""
    rows = parse_yaml(args)
    write_ped(args, rows)


if __name__ == "__main__":
    """ Parse command line arguments """
    parser = argparse.ArgumentParser()

    # Required positional argument
    parser.add_argument("cohortyaml", help="100Humans cohort yaml file", type=str)

    # Required positional argument
    parser.add_argument("cohortid", help="cohort id", type=str)

    # Required positional argument
    parser.add_argument("pedigree", help="output PED file", type=str)

    # Specify output of "--version"
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (version {version})".format(version=__version__))

    args = parser.parse_args()
    main(args)
