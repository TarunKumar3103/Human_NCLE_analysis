#!/usr/bin/env python3
"""Read-only audit of the supplied protein-decay workbook and repository snapshot.

Requires pandas and openpyxl. Does not run notebooks, simulations, or remote jobs.
Exact leading-accession joins only; header-derived alternatives are review flags.
The illustrative pLDDT/length/shape cutoffs reproduce the OLD repository workflow,
not an agreed eligibility protocol for the proposed decay study.
"""
import argparse
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


def fasta_records(path):
    records = {}
    uid = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if line.startswith('>'):
                uid = line.split('|')[1]
                if uid in records:
                    raise ValueError(f'Duplicate FASTA accession: {uid}')
                records[uid] = {'sequence': '', 'header': line}
            elif uid is not None:
                records[uid]['sequence'] += line
    return records


def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--workbook', type=Path, required=True)
    parser.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    d = args.repo / '0_Dataset_curation_and_logistic_regression'
    ent = pd.read_csv(d / 'Human_AF_combined_20250614.csv', sep='|', dtype={'gene': str, 'CCBond': str})
    conf = pd.read_csv(d / 'Human_Avg_pLDDTs.csv').rename(columns={'gene': 'accession', '<pLDDT>': 'mean_plddt', 'std': 'plddt_sd'})
    exp = pd.read_excel(args.workbook, sheet_name='FinalTable_RPE_')
    assert exp['Protein IDs (Uniprot)'].is_unique
    assert conf['accession'].is_unique
    assert set(exp['Degradation profile']) <= {'ED', 'NED', 'UN'}
    bool_map = {'True': True, 'False': False, '1.0': True, '0.0': False}
    unknown = set(ent['CCBond'].dropna().unique()) - set(bool_map)
    if unknown:
        raise ValueError(f'Unexpected CCBond codes: {unknown}')
    ent['covalent'] = ent['CCBond'].map(bool_map)
    assert (ent['ENT-ID'].isna() == ent['CCBond'].isna()).all()
    proteins = []
    for uid, rows in ent.groupby('gene', sort=True):
        calls = rows.loc[rows['ENT-ID'].notna()]
        n_cov = int(calls['covalent'].eq(True).sum())
        n_non = int(calls['covalent'].eq(False).sum())
        status = 'exclude_covalent_lasso' if n_cov else ('NCLE' if n_non else 'no_detected_entanglement')
        sizes = calls['prot_size'].dropna().unique()
        if len(sizes) > 1:
            raise ValueError(f'Conflicting structural lengths: {uid}')
        proteins.append({'accession': uid, 'entanglement_status': status,
                         'n_entanglement_rows': len(calls), 'n_noncovalent_rows': n_non,
                         'n_covalent_rows': n_cov,
                         'structure_length_in_ent_file': sizes[0] if len(sizes) else None,
                         'distinct_ACO_values': int(calls['ACO'].nunique())})
    prot = pd.DataFrame(proteins)
    knots = pd.read_csv(d / 'AF4_knots.dat', sep='\t', dtype=str)
    human_knots = knots.loc[knots['Organism'].eq('Human')].copy()
    human_knots['accession'] = human_knots['Uniprot'].str.split('-').str[0]
    knot_ids = set(human_knots['accession'])
    prot['listed_in_human_knot_file'] = prot['accession'].isin(knot_ids)
    seqs = fasta_records(d / 'UP000005640_9606.fasta')
    membranes = fasta_records(d / 'uniprotkb_human_keywords_membrane_2025_03_19.fasta')
    prot['sequence_length'] = prot['accession'].map({k: len(v['sequence']) for k, v in seqs.items()})
    prot['sequence_length_matches_ent_file'] = pd.Series(pd.NA, index=prot.index, dtype='boolean')
    has_both_lengths = prot['sequence_length'].notna() & prot['structure_length_in_ent_file'].notna()
    prot.loc[has_both_lengths, 'sequence_length_matches_ent_file'] = prot.loc[has_both_lengths, 'sequence_length'].eq(prot.loc[has_both_lengths, 'structure_length_in_ent_file'])
    prot['listed_membrane_keyword_20250319'] = prot['accession'].isin(membranes)
    prot = prot.merge(conf[['accession', 'mean_plddt', 'plddt_sd']], on='accession', how='left', validate='one_to_one')
    shape = pd.read_excel(d / 'AF_Human_shape.xlsx').rename(columns={'Uniprot': 'accession', 'Rg (A)': 'Rg_angstrom', 'Asphericity': 'asphericity'})
    prot = prot.merge(shape, on='accession', how='left', validate='one_to_one')
    exp.insert(0, 'source_excel_row', range(2, len(exp) + 2))
    exp = exp.rename(columns={'Protein IDs (Uniprot)': 'accession'})
    audit = exp.merge(prot, on='accession', how='left', validate='one_to_one', indicator=True)
    audit['exact_leading_accession_match'] = audit.pop('_merge').eq('both')
    audit['header_length'] = audit['Fasta header'].fillna('').str.len()
    audit['header_at_256_character_limit'] = audit['header_length'].eq(256)
    audit['multiple_gene_names'] = audit['Gene names'].fillna('').str.contains(';', regex=False)
    known_ids = set(prot['accession'])
    def candidates(header):
        raw = sorted(set(re.findall(r'(?:sp|tr)\|([^|;\s]+)\|', str(header))))
        bases = sorted({s.split('-')[0] for s in raw} & known_ids)
        return ';'.join(raw), ';'.join(bases), len(bases)
    candidate_data = audit['Fasta header'].map(candidates).tolist()
    for n, name in enumerate(['accessions_in_supplied_header', 'header_base_candidates_for_review', 'n_header_base_candidates']):
        audit[name] = [v[n] for v in candidate_data]
    audit['binary_decay_for_review'] = audit['Degradation profile'].map({'ED': 0, 'NED': 1}).astype('Int64')
    audit['native_entangled_for_review'] = audit['entanglement_status'].map({'NCLE': 1, 'no_detected_entanglement': 0}).astype('Int64')
    flags = []
    for _, row in audit.iterrows():
        reasons = []
        if row['Degradation profile'] == 'UN': reasons.append('undefined_decay')
        if not row['exact_leading_accession_match']:
            reasons.append('no_exact_leading_accession_match')
        else:
            if row['entanglement_status'] == 'exclude_covalent_lasso': reasons.append('covalent_lasso_old_exclusion')
            if row['listed_in_human_knot_file']: reasons.append('knot_list_old_exclusion')
            if row['mean_plddt'] < 85: reasons.append('mean_plddt_below_85_old_cutoff')
            if pd.isna(row['sequence_length']): reasons.append('sequence_missing')
            if row['listed_membrane_keyword_20250319']: reasons.append('membrane_keyword_old_exclusion')
        flags.append(';'.join(reasons))
    audit['illustrative_basic_filter_flags'] = flags
    audit['passes_illustrative_basic_filters'] = audit['illustrative_basic_filter_flags'].eq('')
    audit['additional_review_required'] = 'sequence_and_isoform;per_residue_confidence;disorder;secretory_annotation;global_contact_order;writhe;visual_structure_QC'

    stages = []
    keep = pd.Series(True, index=audit.index)
    filters = [
        ('All experimental entries', keep.copy()),
        ('ED or NED label', audit['Degradation profile'].isin(['ED', 'NED'])),
        ('Exact leading accession in entanglement dataset', audit['exact_leading_accession_match']),
        ('Exclude proteins with any covalent-lasso row', audit['entanglement_status'].isin(['NCLE', 'no_detected_entanglement'])),
        ('Exclude human knot-list entries', ~audit['listed_in_human_knot_file'].eq(True)),
        ('Mean pLDDT >= 85 (old repository cutoff)', audit['mean_plddt'].ge(85)),
        ('Sequence length available', audit['sequence_length'].notna()),
        ('Exclude membrane-keyword list entries', ~audit['listed_membrane_keyword_20250319'].eq(True)),
        ('Length < 800 (old repository cutoff)', audit['sequence_length'].lt(800)),
        ('Asphericity <= 0.1 (old repository cutoff)', audit['asphericity'].le(0.1)),
    ]
    for name, f in filters:
        previous = int(keep.sum())
        keep &= f
        sub = audit.loc[keep]
        counts = sub['Degradation profile'].value_counts().to_dict()
        stages.append({'stage': name, 'n': len(sub), 'removed_this_stage': previous-len(sub),
                       **{c: int(counts.get(c, 0)) for c in ['ED', 'NED', 'UN']},
                       **{f'{e}_{c}': int((sub['entanglement_status'].eq(e) & sub['Degradation profile'].eq(c)).sum())
                          for e in ['NCLE', 'no_detected_entanglement'] for c in ['NED', 'ED']}})
    unmatched = audit.loc[~audit['exact_leading_accession_match']]
    # Missing structures/metadata stay unavailable rather than being imputed as zero.
    length_comparable = prot.dropna(subset=['structure_length_in_ent_file', 'sequence_length'])
    stats = {
        'experimental_shape': [len(exp), exp.shape[1]-1], 'original_experimental_columns': exp.shape[1]-1,
        'profile_counts': exp['Degradation profile'].value_counts().to_dict(),
        'entanglement_rows': len(ent), 'entanglement_unique_proteins': len(prot),
        'protein_entanglement_status_counts': prot['entanglement_status'].value_counts().to_dict(),
        'blank_entanglement_rows': int(ent['ENT-ID'].isna().sum()),
        'confidence_rows': len(conf), 'human_knot_records': len(human_knots), 'human_knot_unique_accessions': len(knot_ids),
        'exact_leading_matches_all_labels': int(audit['exact_leading_accession_match'].sum()),
        'unmatched_leading_accessions': len(unmatched),
        'unmatched_header_candidate_counts': unmatched['n_header_base_candidates'].value_counts().sort_index().to_dict(),
        'headers_at_256_characters': int(audit['header_at_256_character_limit'].sum()),
        'multiple_gene_name_rows': int(audit['multiple_gene_names'].sum()),
        'proteins_with_multiple_ACO_values': int(prot['distinct_ACO_values'].gt(1).sum()),
        'entanglement_vs_fasta_length_disagreements': int(length_comparable['structure_length_in_ent_file'].ne(length_comparable['sequence_length']).sum()),
        'basic_filter_count_not_final_eligibility': int(audit['passes_illustrative_basic_filters'].sum()),
        'audit_policy': 'Only exact leading accessions joined. UN is not a control. Header alternatives are not resolved. Filters are inherited examples, not approved criteria.',
        'source_sha256': {p.name: sha256(p) for p in [args.workbook, d/'Human_AF_combined_20250614.csv', d/'Human_Avg_pLDDTs.csv', d/'AF4_knots.dat', d/'UP000005640_9606.fasta', d/'AF_Human_shape.xlsx', d/'uniprotkb_human_keywords_membrane_2025_03_19.fasta']},
    }
    prot.to_csv(args.out / 'protein_level_structure_metadata.csv', index=False)
    audit.to_csv(args.out / 'human_decay_join_audit.csv', index=False)
    pd.DataFrame(stages).to_csv(args.out / 'illustrative_filter_counts.csv', index=False)
    (args.out / 'audit_summary.json').write_text(json.dumps(stats, indent=2) + '\n')
    print(json.dumps(stats, indent=2))
    print(pd.DataFrame(stages).to_string(index=False))


if __name__ == '__main__':
    main()
