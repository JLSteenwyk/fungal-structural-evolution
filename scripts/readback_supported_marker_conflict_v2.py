#!/usr/bin/env python3
"""Check every diagnostic grid, support status and exhaustive set-based conflict maximum."""
import argparse,json,hashlib
from pathlib import Path
from collections import defaultdict
import pandas as pd
from Bio import Phylo
from audit_joint_path_uncertainty import checked,rows,sha


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    for name in ['diagnostic','snapshot','trees','output']:ap.add_argument('--'+name,type=Path,required=True)
    a=ap.parse_args()
    if a.output.exists():raise FileExistsError(a.output)
    r=checked(a.diagnostic);sr=checked(a.snapshot)
    for name,digest in r['source_pins'].items():
        if sha(Path(name))!=digest:raise ValueError('Changed source')
    support=defaultdict(dict)
    for row in rows(a.snapshot/'branch_support.tsv'):
        support[row['marker']][row['smaller_side_taxa']]=None if row['sh_alrt_percent']=='' else float(row['sh_alrt_percent'])
    data=pd.read_csv(a.diagnostic/'marker_guide_conflict.tsv',sep='\t',keep_default_na=False)
    keys=['marker','guide','full_split_taxa','sh_alrt_cutoff'];assert not data.duplicated(keys).any()
    assert len(data)==r['rows'];assert set(data.marker)=={x['marker'] for x in sr['inputs']}
    sampled=0;all_rows=0
    for (marker,guide),group in data.groupby(['marker','guide']):
        taxa={t.name for t in Phylo.read(a.trees/marker/'tree.treefile','newick').get_terminals()}
        split_sets=[(set(split.split(';')),taxa-set(split.split(';')),score) for split,score in support[marker].items() if score is not None]
        # Independent tree input for the expected full internal-edge grid.
        guide_path=next(Path(x).parent/'guide.treefile' for x in r['source_pins'] if Path(x).parent.name==guide and Path(x).name=='receipt.json')
        gt=Phylo.read(guide_path,'newick');alltaxa={t.name for t in gt.get_terminals()};expected=set()
        for node in gt.find_clades():
            if node is gt.root:continue
            side={t.name for t in node.get_terminals()};parts=[tuple(sorted(side)),tuple(sorted(alltaxa-side))];small=sorted(parts,key=lambda z:(len(z),z))[0]
            if len(small)>1:expected.add(';'.join(small))
        assert set(zip(group.full_split_taxa,group.sh_alrt_cutoff))=={(s,c) for s in expected for c in [80,95]}
        assert group.groupby('full_split_taxa').maximum_conflicting_sh_alrt.nunique(dropna=False).eq(1).all()
        for row in group.to_dict('records'):
            restricted=set(row['full_split_taxa'].split(';'))&taxa
            side=sorted([tuple(sorted(restricted)),tuple(sorted(taxa-restricted))],key=lambda z:(len(z),z))[0]
            assert ';'.join(side)==row['restricted_split_taxa'] and row['marker_taxa']==len(taxa)
            exact=support[marker].get(';'.join(side));observed=row['exact_split_sh_alrt'];observed=float(observed) if observed!='' else '';assert ('' if exact is None else exact)==observed
            maxscore=row['maximum_conflicting_sh_alrt'];maxscore=float(maxscore) if maxscore!='' else '';cutoff=row['sh_alrt_cutoff']
            status='uninformative_taxon_coverage' if len(side)<2 else 'supported_concordance' if exact is not None and exact>=cutoff else 'supported_conflict' if maxscore!='' and float(maxscore)>=cutoff else 'unresolved_at_support_cutoff'
            assert row['status']==status
            if maxscore!='':
                witness=set(row['maximum_conflict_witness_split'].split(';'));x=set(side)
                assert x&witness and x&(taxa-witness) and (taxa-x)&witness and (taxa-x)&(taxa-witness)
                assert support[marker][row['maximum_conflict_witness_split']]==float(maxscore)
            if cutoff==80:
                scores=[];x=set(side)
                if len(x)>=2:
                    nx=taxa-x
                    for y,ny,score in split_sets:
                        if not x.isdisjoint(y) and not x.isdisjoint(ny) and not nx.isdisjoint(y) and not nx.isdisjoint(ny):scores.append(score)
                assert (max(scores) if scores else '')==maxscore
                sampled+=1
            all_rows+=1
        print('Exhaustive set maxima checked',marker,guide,sampled,flush=True)
    assert sampled*2==all_rows
    summary=pd.read_csv(a.diagnostic/'edge_summary.tsv',sep='\t');gkeys=['guide','full_split_taxa','sh_alrt_cutoff'];counts=data.groupby(gkeys+['status']).size().unstack(fill_value=0)
    actual=summary.set_index(gkeys).sort_index();counts=counts.sort_index();assert actual.index.equals(counts.index)
    for status in counts.columns:assert (actual[status]==counts[status]).all()
    assert (actual.completed_markers==sr['completed_markers']).all()
    result={'status':'passed_full_exhaustive_supported_marker_conflict_readback','grid_status_and_witness_rows_checked':all_rows,'independent_full_set_conflict_maxima_checked':sampled,'edge_summary_rows_checked':len(summary),'diagnostic_receipt_sha256':sha(a.diagnostic/'receipt.json'),'script_sha256':sha(Path(__file__)),'scope':'Full guide/marker/cutoff grid, restricted splits, exact supports, classifications, conflict witnesses and aggregate counts checked. Independent exhaustive set-based maximum search on every full-guide edge for every marker/guide; both cutoff classifications share the independently checked maximum.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))

if __name__=='__main__':main()
