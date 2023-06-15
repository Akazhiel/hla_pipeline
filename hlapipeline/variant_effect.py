#!/usr/bin/env python

from Bio.Seq import translate
from Bio.Data import IUPACData
import re


three_to_one = IUPACData.protein_letters_3to1_extended


def translate_dna(seq):
    return translate(seq, to_stop=True)

def retrieve_cDNA(cDNA_dict, effect):
    try:
        cDNA_seq = cDNA_dict[effect.transcript_id]
    except KeyError:
        pass
    return

def retrieve_protein(AA_dict, effect):
    try:
        protein = AA_dict[effect.transcript_id]
    except KeyError:
        protein = effect.original_protein_sequence
    return protein

def missense_variant(starts, ends, wt_mer, mut_mer, errors, mut_dna, mut_aa, transcript, cDNA_pos, aa_pos, effect, cDNA_dict, AA_dict):
    ref_AA, var_AA = [three_to_one[aa] for aa in re.split(r'\d+', mut_aa.lstrip('p.')) if len(aa[1])]
    # cDNA_seq = retrieve_cDNA(cDNA_dict, effect)
    protein_seq = retrieve_protein(AA_dict, effect)
    end = [aa_pos + x if (aa_pos + x) < len(protein_seq) else None for x in ends]
    start = [aa_pos - x for x in starts]
    if any(ele < 0 for ele in start):
        errors += ' Start of sequence is shorter than 12aa from mutation'
        start = [0 if x < 0 else x for x in start]
    wt_mer = [protein_seq[x:y] for x, y in zip(start, end)]
    mut_mer = [protein_seq[x:aa_pos - 1] + var_AA + protein_seq[aa_pos:y] for x, y in zip(start, end)]
    return errors, wt_mer, mut_mer


def inframe_variant(starts, ends, wt_mer, mut_mer, errors, mut_dna, mut_aa, transcript, cDNA_pos, aa_pos, effect, cDNA_dict, AA_dict):
    cDNA_seq = cDNA_dict[transcript]
    protein_seq = AA_dict[transcript]
    end = [aa_pos + x if (aa_pos + x) < len(protein_seq) else None for x in ends]
    start = [aa_pos - x for x in starts]
    if any(ele < 0 for ele in start):
        errors += ' Start of sequence is shorter than 12aa from mutation'
        start = [0 if x < 0 else x for x in start]
    wt_mer = [protein_seq[x:y] for x, y in zip(start, end)]
    if 'dup' in mut_dna:
        dup_pos = list(map(int, re.findall(r'\d+', mut_dna))) if mut_dna != '' else 0
        mut_fasta = cDNA_seq[:dup_pos[0] - 1] + cDNA_seq[dup_pos[0] - 1:dup_pos[-1]] + cDNA_seq[dup_pos[0] - 1:]
        mut_protein = translate_dna(mut_fasta)
        mut_mer = [mut_protein[x:y] for x, y in zip(start, end)]
    elif 'ins' in mut_dna:
        ins_seq = mut_dna[int(mut_dna.find('ins')) + 3:]
        mut_fasta = cDNA_seq[:cDNA_pos] + ins_seq + cDNA_seq[cDNA_pos:]
        mut_protein = translate_dna(mut_fasta)
        mut_mer = [mut_protein[x:y] for x, y in zip(start, end)]
    elif 'del' in mut_dna:
        del_pos = list(map(int, re.findall(r'\d+', mut_dna))) if mut_dna != '' else 0
        mut_fasta = cDNA_seq[:del_pos[0] - 1] + cDNA_seq[del_pos[1]:]
        mut_protein = translate_dna(mut_fasta)
        mut_mer = [mut_protein[x:y] for x, y in zip(start, end)]
    return errors, wt_mer, mut_mer


def frameshift_variant(ref, starts, ends, wt_mer, mut_mer, errors, mut_dna, mut_aa, transcript, cDNA_pos, aa_pos, effect, cDNA_dict, AA_dict, three_prime_utr_dict):
    CDS_seq = cDNA_dict[transcript]
    protein_seq = AA_dict[transcript]
    try:
        utr = three_prime_utr_dict[transcript]
    except KeyError:
        utr = ''
    end = [aa_pos + x if (aa_pos + x) < len(protein_seq) else None for x in ends]
    start = [aa_pos - x for x in starts]
    cDNA_seq = CDS_seq + utr
    if any(ele < 0 for ele in start):
        errors += ' Start of sequence is shorter than 12aa from mutation'
        start = [0 if x < 0 else x for x in start]
    wt_mer = [protein_seq[x:y] for x, y in zip(start, end)]
    if 'del' in mut_dna:
        fs = len(ref)
        mut_cDNA = cDNA_seq[:cDNA_pos - 1] + cDNA_seq[cDNA_pos + fs - 1:]
        mut_fasta = str(translate_dna(mut_cDNA.replace(' ', '')))
        mut_mer = [mut_fasta[x:] for x in start]
    elif 'dup' in mut_dna:
        dup_pos = [None, None]
        dup_pos = list(map(int, re.findall(r'\d+', mut_dna))) if mut_dna != '' else 0
        mut_fasta = cDNA_seq[:dup_pos[0]] + cDNA_seq[dup_pos[0] - 1:dup_pos[-1]] + cDNA_seq[dup_pos[-1]:]
        mut_protein = translate_dna(mut_fasta)
        mut_mer = [mut_protein[x:] for x in start]
    elif 'ins' in mut_dna:
        ins_seq = mut_dna[int(mut_dna.find('ins')) + 3:]
        mut_fasta = cDNA_seq[:cDNA_pos] + ins_seq + cDNA_seq[cDNA_pos:]
        mut_protein = translate_dna(mut_fasta)
        mut_mer = [mut_protein[x:y] for x, y in zip(start, end)]
    return errors, wt_mer, mut_mer


def stoplost_variant(starts, ends, wt_mer, mut_mer, errors, mut_dna, mut_aa, transcript, cDNA_pos, aa_pos, effect, cDNA_dict, AA_dict, three_prime_utr_dict):
    CDS_seq = cDNA_dict[transcript]
    protein_seq = AA_dict[transcript]
    try:
        utr = three_prime_utr_dict[transcript]
    except KeyError:
        utr = ''
    end = [aa_pos + x if (aa_pos + x) < len(protein_seq) else None for x in ends]
    start = [aa_pos - x for x in starts]
    cDNA_seq = CDS_seq + utr
    pass
