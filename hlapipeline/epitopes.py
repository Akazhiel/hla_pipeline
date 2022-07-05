"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""

from varcode import Variant
from Bio.Seq import translate
from Bio.Seq import IUPACData
import re


def translate_dna(seq):
    return translate(seq, to_stop=True)


def create_epitope_varcode(chrm, start, ref, alt, db, mut_dna, mut_aa, transcript, funcensgene, cDNA_dict, AA_dict):
    # Retrieve variant info
    vinfo = Variant(contig=chrm, start=start, ref=ref, alt=alt, ensembl=db, allow_extended_nucleotides=True)
    effect = [effect for effect in vinfo.effects() if effect.transcript_id == transcript][0]
    errors = "Flags:"
    wt_mer = ['-', '-']
    mut_mer = ['-', '-']
    starts = [13, 21]
    ends = [12, 20]
    pos = -1
    three_to_one = IUPACData.protein_letters_3to1_extended
    if effect is None:
        errors += ' could not infer the effect'
    else:
        # Retrieve effect type
        protein_mut = effect.short_description
        if protein_mut is None:
            errors += ' could not retrieve AA mutation'
        elif not protein_mut.startswith('p.'):
            errors += ' invalid mutation {}'.format(protein_mut)
            aa_pos = int(re.findall(r'\d+', mut_aa)[0]) if mut_aa != '' else 0
            cDNA_pos = int(re.findall(r'\d+', mut_dna)[0]) if mut_dna != '' else 0
            if 'missense' in funcensgene:
                ref_AA, var_AA = [three_to_one[aa] for aa in re.split(r'\d+', mut_aa.lstrip('p.'))]
                if aa_pos == 0:
                    errors += ' can not code for this mutated position'
                else:
                    protein_seq = AA_dict[transcript]
                    end = [aa_pos + x if aa_pos + x < len(protein_seq) else None for x in ends]
                    start = [aa_pos - x for x in starts]
                    if any(ele < 0 for ele in start):
                        errors += ' Start of sequence is shorter than 12aa from mutation'
                        start = [0 if x < 0 else x for x in start]
                    wt_mer = [protein_seq[x:y] for x, y in zip(start, end)]
                    mut_mer = [protein_seq[x:aa_pos - 1] + var_AA + protein_seq[aa_pos:y] for x, y in zip(start, end)]
            elif 'inframe' in funcensgene:
                if 'dup' in mut_dna:
                    dup_pos = cDNA_pos - 1
                    dup_base = cDNA_dict[transcript][dup_pos]
                elif 'ins' in mut_dna:
                    pass
                elif 'del' in mut_dna:
                    pass
            elif 'frameshift' in funcensgene:
                fs = len(ref)
                cDNA_seq = cDNA_dict[transcript]
                mut_cDNA = cDNA_seq[:cDNA_pos - 1] + cDNA_seq[cDNA_pos + fs - 1:]
                start = [aa_pos - x for x in starts]
                if any(ele < 0 for ele in start):
                    errors += ' Start of sequence is shorter than 12aa from mutation'
                    start = [0 if x < 0 else x for x in start]
                wt_mer = [effect.original_protein_sequence[x:aa_pos + 13] for x in start]
                mut_fasta = str(translate_dna(mut_cDNA.replace(' ', '')))
                mut_mer = [mut_fasta[x:] for x in start]
        elif protein_mut.startswith('p.X'):
            errors += ' mutation occurs in stop codon'
        else:
            # Retrieve pos
            pos = effect.aa_mutation_start_offset
            protein_seq = effect.original_protein_sequence
            if pos is None:
                errors += ' could not find the position for this mutation'
            elif pos == 0:
                errors += ' mutation occurs in start codon'
            else:
                if effect.mutant_protein_sequence is None or effect.original_protein_sequence is None:
                    errors += ' could not retrieve protein sequence'
                else:
                    # Type of effect
                    effect_type = type(effect).__name__
                    if 'Stop' in effect_type:
                        errors += ' stop mutation'
                    elif 'FrameShift' in effect_type:
                        start = [pos - x + 1 for x in starts]
                        end = [pos + x + 1 if pos + x + 1 < len(protein_seq) else None for x in ends]
                        if any(ele < 0 for ele in start):
                            errors += ' Start of sequence is shorter than 12aa from mutation'
                            start = [0 if x < 0 else x for x in start]
                        wt_mer = [effect.original_protein_sequence[x:y] for x, y in zip(start, end)]
                        mut_mer = [effect.mutant_protein_sequence[x:] for x in start]
                    elif 'Substitution' in effect_type \
                            or 'Deletion' in effect_type:
                        start = [pos - x + 1 for x in starts]
                        end = [pos + x + 1 if pos + x + 1 < len(protein_seq) else None for x in ends]
                        if any(ele < 0 for ele in start):
                            errors += ' Start of sequence is shorter than 12aa from mutation'
                            start = [0 if x < 0 else x for x in start]
                        wt_mer = [effect.original_protein_sequence[x:y] for x, y in zip(start, end)]
                        mut_mer = [effect.mutant_protein_sequence[x:y] for x, y in zip(start, end)]
                    elif 'Insertion' in effect_type:
                        start = [pos - x + 1 for x in starts]
                        if any(ele < 0 for ele in start):
                            errors += ' Start of sequence is shorter than 12aa from mutation'
                            start = [0 if x < 0 else x for x in start]
                        size = int(abs(len(ref) - len(alt)) / 3)
                        wt_mer = [effect.original_protein_sequence[x:pos+size+13] for x in start]
                        mut_mer = [effect.mutant_protein_sequence[x:pos+size+13] for x in start]
                    else:
                        errors += ' unknown exonic function {}'.format(effect_type)
    return pos, errors, wt_mer, mut_mer
