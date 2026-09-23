#!/usr/bin/env python3
"""Account for every unflagged omission after completed native-family replay."""
import argparse
from collections import Counter,defaultdict
import csv
import json
from pathlib import Path
import time
import psutil
from catalog_whole_proteome_structures import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True)
    args=ap.parse_args();plan=json.loads(args.plan.read_text());ph=sha(args.plan)
    def verify():
        if sha(args.plan)!=ph:raise ValueError('Changed plan')
        for p,h in plan['pins'].items():
            if sha(p)!=h:raise ValueError('Changed bound source: '+p)
    verify();out=Path(plan['output'])
    if out.exists():raise FileExistsError(out)
    out.mkdir(parents=True);(out/'state.json').write_text('{"status":"waiting_for_full_replay"}\n')
    while True:
        try:
            p=psutil.Process(plan['predecessor']['pid']);live=p.create_time()==plan['predecessor']['create_time'] and p.status()!=psutil.STATUS_ZOMBIE
        except psutil.NoSuchProcess:live=False
        if not live:break
        time.sleep(20)
    verify();rp=Path(plan['replay'])/'receipt.json';r=json.loads(rp.read_text())
    if r['status']!='passed_all_affected_family_replays' or r['plan_sha256']!=sha(plan['replay_plan']):raise ValueError('Replay incomplete or discrepant')
    data=Path(plan['replay'])/'family_replays.jsonl'
    if sha(data)!=r['artifacts']['family_replays.jsonl']:raise ValueError('Changed replay output')
    expected={};sizes={};flags=Counter()
    with Path(plan['missing']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):
            family=row['family'];sizes[family]=int(row['family_size'])
            if row['native_flagged_misplaced']=='True':flags[family]+=1;continue
            key=(family,row['taxon_id']+'_'+row['protein_id'])
            if key in expected:raise ValueError('Repeated source omission')
            expected[key]=row
    placements={}
    with Path(plan['placements']).open() as f:
        for row in csv.DictReader(f,delimiter='\t'):placements[(row['family'],row['gene_label'])]=row
    if set(placements)!=set(expected):raise ValueError('Independent placement/source grid differs')
    families={k[0] for k in expected};seen=set();genes=set();rows=[];mechanisms=Counter();hogs=0
    with data.open() as f:
        for line in f:
            entry=json.loads(line);family=entry['family']
            if family not in families or family in seen:raise ValueError('Wrong or repeated replay family')
            seen.add(family)
            if entry['status']!='matched' or not all(entry[k] for k in ['root_memberships_match','misplaced_flags_match','omissions_match']):raise ValueError('Unresolved replay mismatch')
            if entry['source_genes']!=sizes[family] or entry['replayed_flagged_genes']!=flags[family]:raise ValueError('Source size/flag count differs')
            positions=entry['omitted_gene_parent_states']
            if len(positions)!=entry['unflagged_missing_genes'] or entry['replayed_missing_genes']!=len(positions)+flags[family]:raise ValueError('Omission counts differ')
            hogs+=entry['root_groups']
            for p in positions:
                key=(family,p['gene'])
                if key not in expected or key in genes:raise ValueError('Wrong/repeated omitted gene')
                genes.add(key);placement=placements[key]
                if p['parent']!=placement['immediate_gene_tree_parent'] or placement['inside_emitted_parent_clade']!='False':raise ValueError('Independent tree position differs')
                classification='root_duplication_blocked_parent_with_leaf_without_group' if p['sp_node']=='N0' and 'N0' in p['dups_below'] else 'other_parent_state_requires_review'
                mechanisms[classification]+=1;source=expected[key]
                rows.append(dict(family=family,taxon_id=source['taxon_id'],protein_id=source['protein_id'],native_gene_id=source['native_gene_id'],parent_node=p['parent'],parent_species_node=p['sp_node'],parent_is_duplication=p['dup'],root_duplication_at_or_below_parent='N0' in p['dups_below'],algorithmic_disposition=classification))
    if seen!=families or genes!=set(expected) or len(seen)!=r['families'] or len(genes)!=r['unflagged_missing_genes']:raise ValueError('Incomplete full replay scope')
    trace=json.loads(Path(plan['placement_receipt']).read_text())
    if hogs!=trace['checked_root_hog_parent_memberships']:raise ValueError('Root group count differs from full placement trace')
    table=out/'omission_dispositions.tsv'
    with table.open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    verify();result={'status':'complete_full_replay_omission_disposition_readback','families':len(seen),'genes':len(genes),'root_hog_memberships_replayed':hogs,'parent_state_dispositions':dict(mechanisms),'replay_receipt_sha256':sha(rp),'plan_sha256':ph,'artifacts':{'omission_dispositions.tsv':sha(table)},'scope':'Every unflagged missing gene accounted for against source IDs and independent tree-position trace; completed native replay reproduced root memberships and misplaced flags. Parent states identify the native writer rule within these resolved trees. Does not independently establish tree rooting, duplication timing, ancestral homology, biological loss or correctness under another guide.'}
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
