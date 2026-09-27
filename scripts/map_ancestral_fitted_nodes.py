#!/usr/bin/env python3
"""Map candidate ancestors by incident tip partitions, never by node label alone."""
import csv,hashlib,json
from collections import Counter,defaultdict
from pathlib import Path
from Bio import Phylo
from prepare_case_ancestral_neighborhoods import sha,read


def signatures(tree):
    allnames={n.name for n in tree.get_terminals()};desc={};result={}
    for node in tree.find_clades(order='postorder'):
        desc[node]=set().union(*(desc[c] for c in node.clades)) if node.clades else {node.name}
        if not node.clades:continue
        parts=[desc[c] for c in node.clades]
        if node is not tree.root:parts.append(allnames-desc[node])
        assert all(parts) and set().union(*parts)==allnames and sum(map(len,parts))==len(allnames)
        result[node]=tuple(sorted(tuple(sorted(p)) for p in parts))
    return result


def main():
    root=Path('results/ancestral/case-local-trees-20260927-v1');tr=root/'receipt.json';r=json.loads(tr.read_text())
    tc=Path('metadata/ancestral_case_local_trees_completed_20260927.json');assert json.loads(tc.read_text())['source_receipt_sha256']==sha(tr)
    nodepath=root/'ancestral_node_mapping.tsv';assert sha(nodepath)==r['artifacts'][nodepath.name]
    nodes=[x for x in read(nodepath) if x['dataset']=='domain']
    fits=Path('results/ancestral/case-domain-model-fits-20260927-v1');fr=fits/'receipt.json';f=json.loads(fr.read_text())
    lr=Path('results/ancestral/case-domain-likelihood-replay-20260927-v1/receipt.json');replay=json.loads(lr.read_text());assert replay['status']=='passed_all_156_independent_domain_likelihood_replays' and replay['source_receipt_sha256']==sha(fr)
    pins={str(p):sha(p) for p in [tr,tc,nodepath,fr,lr]};source={}
    for family in sorted({x['family'] for x in nodes}):
        for guide in ['profile','mafft']:
            p=root/(guide+'-'+family+'-domain.nwk');assert sha(p)==r['artifacts'][p.name];pins[str(p)]=sha(p)
            tree=Phylo.read(p,'newick');sigs=signatures(tree)
            for row in nodes:
                if row['family']!=family or row['guide']!=guide:continue
                names=json.loads(row['output_nodes_json']);assert len(names)==1
                found=[n for n in sigs if n.name==names[0]];assert len(found)==1
                node=found[0];source[family,guide,int(row['level'])]=(row,node.name,sigs[node],node is tree.root)
    rows=[]
    for fit in f['results']:
        job=fit['job'];p=fits/job['job_id']/'fit.treefile';assert sha(p)==fit['artifacts'][p.name];pins[str(p)]=sha(p)
        tree=Phylo.read(p,'newick');sigs=signatures(tree);by=defaultdict(list);parent={c:n for n in tree.find_clades() for c in n.clades}
        for node,sig in sigs.items():by[sig].append(node)
        for guide in ['profile','mafft']:
            for level in range(4):
                original,source_name,sig,isroot=source[job['family'],guide,level]
                matches=by.get(sig,[]);degree=len(sig)
                if degree==2:
                    status='degree_two_root_position_not_identified' if isroot else 'degree_two_internal_position_not_identified'
                    assert not matches
                else:
                    assert degree>=3 and len(matches)==1,(job['job_id'],guide,level,degree,len(matches))
                    status='matched_unrooted_vertex'
                incident=[]
                if matches:
                    node=matches[0];incident=[c.branch_length for c in node.clades]
                    if node in parent:incident.append(node.branch_length)
                    assert len(incident)==degree
                rows.append(dict(job_id=job['job_id'],family=job['family'],boundary=job['boundary'],method=job['method'],model=job['model'],guide=guide,level=level,source_node=source_name,source_is_local_root=int(isroot),incident_partitions=degree,partition_sizes_json=json.dumps([len(s) for s in sig]),partition_signature_sha256=hashlib.sha256(json.dumps(sig,separators=(',',':')).encode()).hexdigest(),status=status,fitted_node=matches[0].name if matches else '',minimum_incident_length=min(incident) if incident else '',incident_lengths_json=json.dumps(incident),retained_descendants=original['retained_descendants']))
    assert len(rows)==1248
    # The two guide aliases must produce the same partitions and disposition,
    # not merely the same node names in a fitted tree.
    grouped=defaultdict(list)
    for row in rows:grouped[row['job_id'],row['level']].append(row)
    for pair in grouped.values():
        assert len(pair)==2
        for key in ['partition_signature_sha256','status','fitted_node']:assert pair[0][key]==pair[1][key]
    out=Path('results/ancestral/case-fitted-node-mapping-20260927-v1');out.mkdir(exist_ok=False)
    with (out/'node_mapping.tsv').open('w') as handle:w=csv.DictWriter(handle,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
    assert read(out/'node_mapping.tsv')==[{k:str(v) for k,v in row.items()} for row in rows]
    for p,h in pins.items():assert sha(p)==h
    result=dict(status='complete_candidate_ancestor_mapping_by_incident_tip_partitions',fits=156,guide_node_rows=len(rows),unique_fit_node_candidates=len(grouped),status_counts=dict(Counter(x['status'] for x in rows)),source_hashes=pins,script_sha256=sha(__file__),artifacts={'node_mapping.tsv':sha(out/'node_mapping.tsv')},scope='Exact full incident tip partitions identify unrooted vertices independently of labels across all fits and both guide aliases. Degree-two source-root positions are explicitly unidentified after unrooted fitting. Mapping is not evidence of ancestral polarity, model convergence, node support or a posterior sequence.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='source_hashes'},indent=2))


if __name__=='__main__':main()
