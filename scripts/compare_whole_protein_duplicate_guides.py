#!/usr/bin/env python3
"""Compare exact duplicate gene pairs across guides, retaining missingness and all screens."""
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np
import pandas as pd
from screen_duplication_domain_alignment_coverage import sha


def main():
    root = Path('results/structural_comparisons/whole-protein-reference-sensitivity-20260927-v2')
    rp = root/'receipt.json'; receipt = json.loads(rp.read_text())
    assert receipt['status'] == 'complete_descriptive_whole_protein_reference_sensitivity'
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest
    p = root/'event_sensitivity.tsv'
    d = pd.read_csv(p, sep='\t', keep_default_na=False)
    keys = ['family','gene_a','gene_b','screen']
    assert set(d.guide) == {'mafft','profile'} and len(d) == 209454
    # The upstream producer fixes lexical A/B orientation; fail rather than silently flip it.
    assert (d.gene_a < d.gene_b).all() and not d.duplicated(['guide']+keys).any()
    halves = {g: d[d.guide.eq(g)].drop(columns='guide') for g in ['mafft','profile']}
    combined = halves['mafft'].merge(halves['profile'], on=keys, how='outer',
        suffixes=('_mafft','_profile'), validate='one_to_one', indicator='guide_membership')
    combined['guide_membership'] = combined.guide_membership.astype(str)
    combined = combined.fillna('')
    complete = combined.guide_membership.eq('both')
    eligible = complete & combined.all_references_eligible_mafft.eq(1) & combined.all_references_eligible_profile.eq(1)
    sign = combined.structural_direction_mafft
    stable = eligible & sign.isin(['a_farther','b_farther']) & sign.eq(combined.structural_direction_profile)
    combined['both_guides_all_references_eligible'] = eligible.astype(int)
    combined['same_structural_direction_both_guides'] = stable.astype(int)
    combined['sequence_direction_agrees_both_guides'] = (eligible & combined.sequence_direction_mafft.eq(combined.sequence_direction_profile)).astype(int)
    combined['sequence_structure_same_direction_both_guides'] = (stable & sign.eq(combined.sequence_direction_mafft) & sign.eq(combined.sequence_direction_profile)).astype(int)
    # Full independent dictionary join validates every retained source cell and membership.
    dictionaries = {g:{} for g in ['mafft','profile']}
    with p.open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            key = tuple(row[k] for k in keys); guide = row.pop('guide')
            assert key not in dictionaries[guide]; dictionaries[guide][key] = row
    assert len(combined) == len(set(dictionaries['mafft']) | set(dictionaries['profile']))
    counts = Counter()
    for row in combined.to_dict('records'):
        key = tuple(row[k] for k in keys)
        membership = 'both' if all(key in dictionaries[g] for g in dictionaries) else 'left_only' if key in dictionaries['mafft'] else 'right_only'
        assert row['guide_membership'] == membership
        for guide in dictionaries:
            source = dictionaries[guide].get(key)
            for column in halves[guide].columns:
                if column in keys: continue
                actual = row[column+'_'+guide]
                wanted = source[column] if source else ''
                if wanted == '': assert actual == ''
                elif column.endswith(('_min','_max')): assert np.isclose(float(actual),float(wanted),rtol=1e-12,atol=1e-12)
                elif column in ['reference_links','eligible_reference_links','all_references_eligible']: assert float(actual)==int(wanted)
                else: assert str(actual)==wanted
        a = dictionaries['mafft'].get(key); b = dictionaries['profile'].get(key)
        all_ok = bool(a and b and a['all_references_eligible']=='1' and b['all_references_eligible']=='1')
        same = bool(all_ok and a['structural_direction'] in ['a_farther','b_farther'] and a['structural_direction']==b['structural_direction'])
        assert row['both_guides_all_references_eligible']==all_ok and row['same_structural_direction_both_guides']==same
        counts[(row['screen'],membership)]+=1
    summaries = []
    for screen, group in combined.groupby('screen',sort=True):
        robust = group[group.same_structural_direction_both_guides.eq(1)]
        assert (robust.gene_a.str.split('_').str[0] == robust.gene_b.str.split('_').str[0]).all()
        summaries.append(dict(screen=screen, pair_union=len(group),pairs_in_both=int(group.guide_membership.eq('both').sum()),
            mafft_only=counts[screen,'left_only'],profile_only=counts[screen,'right_only'],
            all_references_eligible_both=int(group.both_guides_all_references_eligible.sum()),
            same_structural_direction_both=len(robust),
            same_sequence_and_structural_direction_both=int(group.sequence_structure_same_direction_both_guides.sum()),
            robust_families=robust.family.nunique(),robust_taxa=robust.gene_a.str.split('_').str[0].nunique()))
    out = Path('results/structural_comparisons/whole-protein-cross-guide-sensitivity-20260927-v1');out.mkdir(exist_ok=False)
    artifacts = {}
    for name,data in [('all_pair_comparisons.tsv',combined),('summary.tsv',pd.DataFrame(summaries))]:
        path=out/name;data.to_csv(path,sep='\t',index=False)
        back=pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False)
        # Explicit string serialization verifies all columns, including missing guide values.
        assert back.equals(data.astype(str).reset_index(drop=True))
        artifacts[name]=sha(path)
    assert sha(p)==receipt['artifacts'][p.name]
    result=dict(status='complete_exact_gene_pair_cross_guide_sensitivity',source_hashes={str(rp):sha(rp),str(p):sha(p)},
        script_sha256=sha(__file__),source_event_screen_rows=len(d),cross_guide_pair_screen_rows=len(combined),summaries=summaries,artifacts=artifacts,
        scope='Exact family and lexical gene-pair matching; guide-specific nodes and missing pairs retained. Independent full dictionary outer join checked. Same gene pairs across guides are not independent events. Sign robustness uses numerical tolerance and is not significant structural asymmetry, selection or excess structure given sequence.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__': main()
