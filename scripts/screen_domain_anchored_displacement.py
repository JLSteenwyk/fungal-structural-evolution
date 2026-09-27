#!/usr/bin/env python3
"""Reapply exact original-interval coverage screens to domain anchors and outside sets."""
import csv,json
from collections import defaultdict,Counter
from fractions import Fraction
from pathlib import Path
from screen_duplication_domain_alignment_coverage import sha


def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    base=Path('results/structural_comparisons');sources={}
    def checked(root,name):
        rp=root/'receipt.json';r=json.loads(rp.read_text());p=root/name
        assert sha(p)==r['artifacts'][name];sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    root=base/'domain-anchored-displacement-20260927-v1'
    measurements=rows(checked(root,'pair_displacements.tsv'))
    partitions=[json.loads(x) for x in checked(root,'residue_partitions.jsonl').open()]
    maps=base/'duplication-domain-common-residues-20260927-v1'
    triads={r['triad_key']:r for r in map(json.loads,checked(maps,'triads.jsonl').open())}
    inputs={}
    for line in checked(base/'duplication-domain-inputs-20260926-v1','inputs.jsonl').open():
        r=json.loads(line)
        if r['mask']=='full':inputs[r['interval_id']]=r
    fields=['domain_triad','whole_triad','mask','order_ab','order_ar','order_br','mapping_definition']
    key=lambda r:tuple(str(r[k]) for k in fields)
    pairgroups=defaultdict(dict)
    for r in measurements:
        assert r['pair'] not in pairgroups[key(r)];pairgroups[key(r)][r['pair']]=r
    screens=[(n,c) for n in [30,50] for c in [50,70,90]];screened=[]
    for part in partitions:
        k=key(part);group=pairgroups[k];assert set(group)=={'ab','ar','br'}
        ds=triads[part['domain_triad']];bb=[inputs[ds['interval_'+role]] for role in ['a','b','reference']]
        lengths=[b['end']-b['start']+1 for b in bb];outside_lengths=[b['original_length']-l for b,l in zip(bb,lengths)]
        ni=len(part['inside_triples']);no=len(part['outside_triples'])
        assert all(int(r['inside_residues'])==ni and int(r['outside_residues'])==no for r in group.values())
        assert all(all(b['start']<=t[i]<=b['end'] for i,b in enumerate(bb)) for t in part['inside_triples'])
        assert all(all(not b['start']<=t[i]<=b['end'] for i,b in enumerate(bb)) for t in part['outside_triples'])
        for n,c in screens:
            anchor_ok=ni>=n and all(100*ni>=c*l for l in lengths) and all(r['status']=='computed' for r in group.values())
            outside_ok=no>=n and all(l>0 and 100*no>=c*l for l in outside_lengths)
            # Rational arithmetic independently checks threshold inclusion and denominators.
            assert anchor_ok==(ni>=n and min(Fraction(ni,l) for l in lengths)>=Fraction(c,100) and all(r['status']=='computed' for r in group.values()))
            assert outside_ok==(no>=n and all(l>0 and Fraction(no,l)>=Fraction(c,100) for l in outside_lengths))
            why=[]
            if ni<n:why.append('short_anchor')
            if any(100*ni<c*l for l in lengths):why.append('low_original_interval_coverage')
            if any(r['status']!='computed' for r in group.values()):why.append('unusable_coordinate_measurement')
            screened.append(dict(zip(fields,k),screen=f'n{n}_c{c}',anchor_residues=ni,outside_residues=no,
                interval_lengths_json=json.dumps(lengths),original_outside_lengths_json=json.dumps(outside_lengths),
                anchor_pass=int(anchor_ok),outside_pass=int(outside_ok),anchor_and_outside_pass=int(anchor_ok and outside_ok),anchor_exclusions=';'.join(why)))
    assert len(screened)==832*6
    case=base/'whole-domain-case-dossiers-20260927-v1';dossiers=rows(checked(case,'case_dossiers.tsv'));links=rows(checked(case,'domain_reference_links.tsv'))
    casekeys=['family','gene_a','gene_b','pfam_accession'];cases=defaultdict(set)
    for r in links:cases[tuple(r[k] for k in casekeys)].add(r['triad_key'])
    index=defaultdict(list)
    for r in screened:index[r['domain_triad'],r['screen']].append(r)
    output=[]
    for case_row in dossiers:
        ck=tuple(case_row[k] for k in casekeys);tks=cases[ck];assert tks
        for n,c in screens:
            screen=f'n{n}_c{c}';pool=[r for tk in tks for r in index[tk,screen]]
            assert len(pool)==32*len(tks)
            output.append(dict(zip(casekeys,ck),species_name=case_row['species_name'],screen=screen,
                alternatives=len(pool),anchor_passing=sum(r['anchor_pass'] for r in pool),anchor_and_outside_passing=sum(r['anchor_and_outside_pass'] for r in pool),
                all_anchors_pass=int(all(r['anchor_pass'] for r in pool)),all_anchors_and_outside_pass=int(all(r['anchor_and_outside_pass'] for r in pool)),
                minimum_anchor_residues=min(r['anchor_residues'] for r in pool),minimum_outside_residues=min(r['outside_residues'] for r in pool)))
    assert len(output)==78
    summary=[]
    for n,c in screens:
        screen=f'n{n}_c{c}';pool=[r for r in output if r['screen']==screen]
        summary.append(dict(screen=screen,cases=len(pool),all_anchors_pass=sum(r['all_anchors_pass'] for r in pool),all_anchors_and_outside_pass=sum(r['all_anchors_and_outside_pass'] for r in pool)))
    out=base/'domain-anchored-coverage-20260927-v1';out.mkdir(exist_ok=False);artifacts={}
    for name,data in [('partition_screens.tsv',screened),('case_screens.tsv',output),('summary.tsv',summary)]:
        p=out/name
        with p.open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t');w.writeheader();w.writerows(data)
        assert rows(p)==[{k:str(v) for k,v in r.items()} for r in data];artifacts[name]=sha(p)
    for path,digest in sources.items():assert sha(path)==digest
    receipt=dict(status='complete_domain_anchor_coverage_screens',source_hashes=sources,script_sha256=sha(__file__),partition_screen_rows=len(screened),case_screen_rows=len(output),summaries=summary,artifacts=artifacts,
        scope='Six exact length/coverage screens applied to all three original domain intervals; outside coverage separately uses full original protein length minus interval length. Both masks, mappings, native orders and interval definitions required for all-alternative flags. Passing is descriptive eligibility, not uncertainty calibration or biological significance.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
