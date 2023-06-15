#!/usr/bin/env python

from varcode import Variant
from hlapipeline.variant_effect import *
import re


def create_epitope_varcode(chrm, start, ref, alt, db, mut_dna, mut_aa, transcript, funcensgene, cDNA_dict, AA_dict, three_prime_utr_dict):
    # Retrieve variant info
    vinfo = Variant(contig=chrm, start=start, ref=ref, alt=alt, ensembl=db, allow_extended_nucleotides=True)
    effect = vinfo.effect_on_transcript(db.transcript_by_id(transcript))
    errors = "Flags:"
    wt_mer = ['-', '-']
    mut_mer = ['-', '-']
    starts = [13, 21]
    ends = [12, 20]
    pos = -1
    if effect is None:
        errors += ' Could not infer the effect.'
    else:
        # Retrieve effect type
        protein_mut = effect.short_description
        if protein_mut is None:
            errors += ' Could not retrieve AA mutation.'
        elif not protein_mut.startswith('p.'):
            errors += ' Computed with dictionary method.'
            errors += ' Invalid mutation {}.'.format(protein_mut)
            aa_pos = int(re.findall(r'\d+', mut_aa)[0]) if mut_aa != '' else 0
            cDNA_pos = int(re.findall(r'\d+', mut_dna)[0]) if mut_dna != '' else 0
            cDNA_pos = cDNA_pos + 1 if cDNA_pos == 1 else cDNA_pos
            if aa_pos == 0:
                errors += ' Can not code for this mutated position.'
            else:
                if 'missense' in funcensgene:
                    errors, wt_mer, mut_mer = missense_variant(starts, ends, wt_mer, mut_mer, errors,
                                                                    mut_dna, mut_aa, transcript, cDNA_pos,
                                                                    aa_pos, effect, cDNA_dict, AA_dict)
                elif 'frameshift' in funcensgene:
                    errors, wt_mer, mut_mer = frameshift_variant(ref, starts, ends, wt_mer, mut_mer, errors,
                                                                    mut_dna, mut_aa, transcript, cDNA_pos,
                                                                    aa_pos, effect, cDNA_dict, AA_dict, three_prime_utr_dict)
                elif 'stop_lost' in funcensgene:
                    errors += ' Stop lost not yet implemented.'
                    pass
                    # errors, wt_mer, mut_mer = stoplost_variant(starts, ends, wt_mer, mut_mer, errors,
                    #                                                 mut_dna, mut_aa, transcript, cDNA_pos,
                    #                                                 aa_pos, effect, cDNA_dict, AA_dict, three_prime_utr_dict)
                elif 'inframe' in funcensgene:
                    errors, wt_mer, mut_mer = inframe_variant(starts, ends, wt_mer, mut_mer, errors,
                                                                mut_dna, mut_aa, transcript, cDNA_pos,
                                                               aa_pos, effect, cDNA_dict, AA_dict)
        elif protein_mut.startswith('p.X'):
            errors += ' Mutation occurs in stop codon.'
        else:
            # Retrieve pos
            pos = effect.aa_mutation_start_offset
            protein_seq = effect.original_protein_sequence
            if pos is None:
                errors += ' Could not find the position for this mutation.'
            elif pos == 0:
                errors += ' Mutation occurs in start codon.'
            else:
                if effect.mutant_protein_sequence is None or effect.original_protein_sequence is None:
                    errors += ' Could not retrieve protein sequence.'
                else:
                    errors += ' Computed with varcode method.'
                    # Type of effect
                    effect_type = type(effect).__name__
                    if 'StopLoss' in effect_type:
                        start = [pos - x + 1 for x in starts]
                        end = [pos + x + 1 if pos + x + 1 < len(protein_seq) else None for x in ends]
                        if any(ele < 0 for ele in start):
                            errors += ' Start of sequence is shorter than 12aa from mutation.'
                            start = [0 if x < 0 else x for x in start]
                        wt_mer = [effect.original_protein_sequence[x:y] for x, y in zip(start, end)]
                        mut_mer = [effect.mutant_protein_sequence[x:] for x in start]
                    elif 'Stop' in effect_type:
                        errors += ' Stop mutation.'
                    elif 'FrameShift' in effect_type:
                        start = [pos - x + 1 for x in starts]
                        end = [pos + x + 1 if pos + x + 1 < len(protein_seq) else None for x in ends]
                        if any(ele < 0 for ele in start):
                            errors += ' Start of sequence is shorter than 12aa from mutation.'
                            start = [0 if x < 0 else x for x in start]
                        wt_mer = [effect.original_protein_sequence[x:y] for x, y in zip(start, end)]
                        mut_mer = [effect.mutant_protein_sequence[x:] for x in start]
                    elif 'Substitution' in effect_type \
                            or 'Deletion' in effect_type:
                        start = [pos - x + 1 for x in starts]
                        end = [pos + x + 1 if pos + x + 1 < len(protein_seq) else None for x in ends]
                        if any(ele < 0 for ele in start):
                            errors += ' Start of sequence is shorter than 12aa from mutation.'
                            start = [0 if x < 0 else x for x in start]
                        wt_mer = [effect.original_protein_sequence[x:y] for x, y in zip(start, end)]
                        mut_mer = [effect.mutant_protein_sequence[x:y] for x, y in zip(start, end)]
                    elif 'Insertion' in effect_type:
                        start = [pos - x + 1 for x in starts]
                        if any(ele < 0 for ele in start):
                            errors += ' Start of sequence is shorter than 12aa from mutation.'
                            start = [0 if x < 0 else x for x in start]
                        size = int(abs(len(ref) - len(alt)) / 3)
                        wt_mer = [effect.original_protein_sequence[x:pos+size+13] for x in start]
                        mut_mer = [effect.mutant_protein_sequence[x:pos+size+13] for x in start]
                    else:
                        errors += ' Unknown exonic function {}.'.format(effect_type)
    return errors, wt_mer, mut_mer
