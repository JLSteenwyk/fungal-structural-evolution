#!/usr/bin/env python3
"""Compare every validated ML split among the first three crossed PMSF runs."""
import csv
import hashlib
import itertools
import json
from pathlib import Path
from Bio import Phylo

RUNS = {
 'profile_profile': ('pmsf-profile-profile-recovery-v2','pmsf-profile-profile-readback-v1'),
 'profile_mafft': ('pmsf-profile-mafft-v1','pmsf-profile-mafft-readback-v1'),
 'mafft_profile': ('pmsf-mafft-profile-v1','pmsf-mafft-profile-readback-20260927-v1')}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main():
    pins={};sources={};universe=None
    for label,(run_name,audit_name) in RUNS.items():
        run=Path('results/phylogeny')/run_name;audit=Path('results/phylogeny')/audit_name
        proof=json.loads((audit/'receipt.json').read_text());source=json.loads((run/'receipt.json').read_text())
        assert proof['status']=='passed_pmsf_profile_tree_and_bootstrap_readback'
        assert proof['source_receipt_sha256']==sha(run/'receipt.json')
        for root,receipt in [(run,source),(audit,proof)]:
            pins[str(root/'receipt.json')]=sha(root/'receipt.json')
            for name,digest in receipt['artifacts'].items():assert sha(root/name)==digest
        taxa={n.name for n in Phylo.read(run/'pmsf.treefile','newick').get_terminals()}
        assert len(taxa)==526
        if universe is None:universe=taxa
        assert taxa==universe
        with (audit/'branch_support.tsv').open() as f:rows=[r for r in csv.DictReader(f,delimiter='\t') if r['tree']=='ml']
        splits={}
        for r in rows:
            side=set(json.loads(r['split_taxa_json']))
            key=min(tuple(sorted(side)),tuple(sorted(taxa-side)),key=lambda x:(len(x),x))
            assert key not in splits and 1<len(key)<525
            splits[key]=r
        assert len(splits)==523
        sources[label]=splits
    comparison=[];conflicts=[];union=[]
    for split in sorted(set().union(*(set(s) for s in sources.values()))):
        for label,splits in sources.items():
            r=splits.get(split)
            union.append(dict(run=label,split_taxa_json=json.dumps(split),present=r is not None,branch_length=r['branch_length'] if r else '',sh_alrt=r['sh_alrt_percent'] if r else '',empirical_ufb=r['empirical_ufboot_percent'] if r else ''))
    for a,b in itertools.combinations(sources,2):
        left,right=sources[a],sources[b];common=set(left)&set(right);nconf=nhigh=0
        for x,y in itertools.product(set(left)-common,set(right)-common):
            xs,ys=set(x),set(y)
            cells=[xs&ys,xs-ys,ys-xs,universe-(xs|ys)]
            if not all(cells):continue
            nconf+=1
            av,bv=left[x],right[y]
            high=all(float(r['sh_alrt_percent'])>=80 and float(r['empirical_ufboot_percent'])>=95 for r in (av,bv))
            nhigh+=high
            conflicts.append(dict(run_a=a,run_b=b,split_a_taxa_json=json.dumps(x),split_b_taxa_json=json.dumps(y),sh_alrt_a=av['sh_alrt_percent'],sh_alrt_b=bv['sh_alrt_percent'],empirical_ufb_a=av['empirical_ufboot_percent'],empirical_ufb_b=bv['empirical_ufboot_percent'],both_meet_existing_support_label=high,witness_quartet=';'.join(min(c) for c in cells)))
        comparison.append(dict(run_a=a,run_b=b,shared_internal_splits=len(common),unique_internal_splits_a=len(left)-len(common),unique_internal_splits_b=len(right)-len(common),rf_distance=len(left)+len(right)-2*len(common),normalized_rf=(len(left)+len(right)-2*len(common))/1046,incompatible_split_pairs=nconf,both_supported_incompatible_pairs=nhigh))
    root=Path('results/phylogeny/pmsf-three-run-topology-sensitivity-20260927-v1');root.mkdir(parents=True,exist_ok=False)
    for name,rows in [('comparisons.tsv',comparison),('split_presence.tsv',union),('conflicts.tsv',conflicts)]:
        with (root/name).open('w') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    for path,digest in pins.items():assert sha(path)==digest
    receipt=dict(status='complete_three_run_pmsf_topology_comparison_pending_readback',taxa=526,runs=RUNS,comparisons=comparison,shared_all_three=len(set.intersection(*(set(s) for s in sources.values()))),sources=pins,script_sha256=sha(__file__),artifacts={p.name:sha(p) for p in root.iterdir()},scope='All ML splits and incompatible unique-split pairs retained. Existing high-support label means SH-aLRT >=80 and empirical UFBoot >=95 on both splits; not correctness or significance. Different alignments have different site universes; neither raw likelihood nor site-profile comparisons are used. Fourth crossed run, consensus sensitivity, rooting, adequacy and discordance remain pending.')
    (root/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps({k:receipt[k] for k in ['comparisons','shared_all_three']},indent=2))


if __name__=='__main__':main()
