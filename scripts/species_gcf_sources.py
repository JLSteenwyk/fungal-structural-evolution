"""Full audited 250-marker/four-supported-baseline source I/O for native gCF."""
import json
from pathlib import Path
from Bio import Phylo
from four_run_pmsf_sources import load_sources
from reference_measurement_union_sources import bind,verify
from run_ortholog_pair_guide_comparison import sha


def load(plan,plan_path):
    bindings=dict(plan['pins']);bind(bindings,plan_path)
    source_plan=json.loads(Path(plan['species_plan']).read_text())
    baseline=json.loads(Path(plan['species_completion']).read_text())
    assert baseline['status']=='complete_verified_full_four_run_pmsf_ML_consensus_sensitivity' and len(baseline['services'])==2
    assert baseline['source_hashes'][plan['species_plan']]==sha(plan['species_plan'])
    sources,manifest,prior=load_sources(source_plan,plan['species_plan'])
    for record in [baseline,dict(source_hashes=prior)]:
        for p,d in record['source_hashes'].items():bind(bindings,p,d)
    bind(bindings,plan['species_completion'])
    marker_closed=json.loads(Path(plan['marker_completion']).read_text())
    assert marker_closed['status']=='complete_verified_full_mafft_marker_support_and_alignment_topology' and len(marker_closed['services'])==3
    assert marker_closed['summary']['markers']==125
    for p,d in marker_closed['source_hashes'].items():bind(bindings,p,d)
    bind(bindings,plan['marker_completion'])
    universe={r['taxon_id'] for r in manifest};markers={}
    for label,spec in plan['markers'].items():
        srp=Path(spec['support'])/'receipt.json';arp=Path(spec['audit'])/'receipt.json'
        sr,ar=[json.loads(p.read_text()) for p in [srp,arp]]
        assert sr['status']=='complete_support_audit' and sr['completed_markers']==sr['planned_markers']==125 and not sr['pending_markers']
        assert ar['status']=='passed_all_snapshot_input_and_graph_split_readbacks' and ar['markers']==125 and ar['snapshot_receipt_sha256']==sha(srp)
        for root,record in [(Path(spec['support']),sr),(Path(spec['audit']),ar)]:
            bind(bindings,root/'receipt.json')
            for name,digest in record['artifacts'].items():bind(bindings,root/name,digest)
        items={r['marker']:r for r in sr['inputs']};assert len(items)==len(sr['inputs'])==125;collection=[]
        for marker,proof in sorted(items.items()):
            root=Path(spec['trees'])/marker;rp=root/'receipt.json';tp=root/'tree.treefile';r=json.loads(rp.read_text())
            bind(bindings,rp,proof['receipt_sha256']);bind(bindings,tp,proof['tree_sha256'])
            assert r['returncode']==0 and r['tree_sha256']==sha(tp)
            command=r['command'];input_path=command[command.index('-s')+1];bind(bindings,input_path,r['input_sha256'])
            tree=Phylo.read(tp,'newick');tips=[t.name for t in tree.get_terminals()]
            assert len(tips)==len(set(tips))==r['taxa'] and set(tips)<=universe and len(tips)>=4
            assert len(list(tree.find_clades()))-1==2*len(tips)-3
            collection.append(dict(marker=marker,path=str(tp),tree_sha256=sha(tp),taxa=len(tips)))
        markers[label]=collection
    assert set(markers)=={'profile','mafft'} and [r['marker'] for r in markers['profile']]==[r['marker'] for r in markers['mafft']]
    verify(bindings)
    return sources,universe,markers,bindings
