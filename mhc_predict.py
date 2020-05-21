#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: jfnavarro
"""
from hlapipeline.common import exec_command
from collections import Counter
from _collections import defaultdict
import json
import argparse


def compute_MHC(hla_exome_cancer, hla_exome_normal, hla_rna, overlap_final):
    HLA_dict = defaultdict(list)

    # First parse hla_exome_cancer and normal (HLA-LA format)
    print('Loading DNA tumor HLAs..')
    for file in hla_exome_cancer:
        with open(file) as f:
            for line in f.readlines():
                columns = line.strip().split('\t')
                hla = columns[1].split("_")[-1]
                alleles = columns[2:]
                HLA_dict[hla].extend(alleles)

    print('Loading DNA normal HLAs..')
    for file in hla_exome_normal:
        with open(file) as f:
            for line in f.readlines():
                columns = line.strip().split('\t')
                hla = columns[1].split("_")[-1]
                alleles = columns[2:]
                HLA_dict[hla].extend(alleles)

    # Parse RNA hlas (arcasHLA JSON format)
    print('Loading RNA HLAs..')
    for file in hla_rna:
        with open(file) as f:
            local_dict = json.load(f)
            for hla,alleles in local_dict.items():
                HLA_dict[hla].extend(alleles)

    # Filter HLAs by occurrences
    filtered_hla = []
    for alleles in HLA_dict.values():
        f = Counter(alleles)
        items = f.most_common()
        cutoff = items[0][1] if len(items) == 1 else items[1][1]
        filtered_hla += [y[0] for y in items if y[1] >= cutoff]
    # MHCflurry format (HLA-A*02:01)
    filtered_hla = ['HLA-{}'.format(':'.join(x.split(':')[0:2])) for x in filtered_hla]
    print('Alleles found: {}'.format(','.join(filtered_hla)))

    # Filter HLAs by allowed alleles in MHCflurry
    allowed_alleles = set()
    with open("alleles.txt") as f:
        for x in f.readlines():
            allowed_alleles.add(x.strip())
    filtered_hla = [x for x in filtered_hla if x in allowed_alleles]
    print('Alleles allowed in MHCflurry: {}'.format(','.join(filtered_hla)))

    # Create protein FASTA file
    print('Creating protein sequencess..')
    added_proteins = set()
    with open('protein_sequences.fasta', 'w') as fwrite:
        with open(overlap_final, 'r') as fread:
            lines = fread.readlines()
            header = lines.pop(0).strip().split('\t')
            for line in lines:
                columns = line.strip().split('\t')
                if int(columns[header.index('Number of Exomes samples (passing)')]) > 0:
                    protein_name = '{} {} {} {}'.format(columns[header.index('Variant key')],
                                                        columns[header.index('transcript ID')],
                                                        columns[header.index('cDNA change')],
                                                        columns[header.index('AA change')])
                    protein_seq = columns[header.index('Mut Epitope')].strip().replace('*','')
                    if protein_seq != '-' and protein_seq not in added_proteins:
                        fwrite.write('>{}\n{}\n'.format(protein_name, protein_seq))
                        added_proteins.add(protein_seq)


    # Run prediction
    print('Predicting MHCs..')
    cmd = 'mhcflurry-predict-scan protein_sequences.fasta --alleles {} ' \
          '--results-all --out predictions.csv --peptide-lengths 8 9 10 11 12'.format(' '.join(filtered_hla))
    exec_command(cmd)
    print('Completed')

parser = argparse.ArgumentParser(description='Script to predict MHCs using the data from the exome and rnaseq variant calling pipelines '
                                             '(created by Jose Fernandez <jc.fernandes.navarro@gmail.com>)',
                                 prog='mhc_predict.py',
                                 usage='mhc_predict.py [options] '
                                       '--hla-dna-normal [file/s with HLA predictions from DNA (Normal)] '
                                       '--hla-dna-tumor [file/s with HLA predictions from DNA (Tumor)] '
                                       '--hla-rna [file/s with HLA predictions from RNA]'
                                       '--variants [file with the final variants generated with merge_results.py]')

parser.add_argument('--hla-dna-normal', nargs='+', default=None, required=True,
                    help='A file or files containing predicted HLAs from normal DNA (table format)')
parser.add_argument('--hla-dna-tumor', nargs='+', default=None, required=True,
                    help='A file or files containing predicted HLAs from tumor DNA (table format)')
parser.add_argument('--hla-rna', nargs='+', default=None, required=True,
                    help='A file or files containing predicted HLAs from RNA (JSON format)')
parser.add_argument('--variants', default=None, required=True,
                    help='A file with the final variants generated with merge_results.py (table format)')

args = parser.parse_args()
compute_MHC(args.hla_dna_tumor, args.hla_dna_normal, args.hla_rna, args.variants)