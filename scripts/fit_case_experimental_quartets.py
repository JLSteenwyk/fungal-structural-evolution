#!/usr/bin/env python3
"""Fit all sequence-anchored experimental quartets on identical observed residues."""
import csv,gzip,json,subprocess,time,hashlib
from collections import defaultdict
from pathlib import Path
import numpy as np
import psutil
from duplication_alignment_numeric_readback import load_pdb
from assess_domain_alignment_geometry import geometry
from measure_domain_anchored_displacement import transform,residual,fixture
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results')
OUTPUT=BASE/'experimental_structures/whole-domain-case-quartet-fits-20260927-v1'
ROLES=['a','b','reference']
PAIRS=[('ab',0,1),('ar',0,2),('br',1,2),('ae',0,3),('be',1,3),('re',2,3)]


def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))


def main():
    fixtures=fixture();sources={}
    for name in ['duplication_alignment_numeric_readback.py','assess_domain_alignment_geometry.py','measure_domain_anchored_displacement.py']:
        p=Path('scripts')/name;sources[str(p)]=sha(p)
    def checked(root,name):
        rp=root/'receipt.json';p=root/name;r=json.loads(rp.read_text());assert sha(p)==r['artifacts'][name]
        sources[str(rp)]=sha(rp);sources[str(p)]=sha(p);return p
    lp=Path('metadata/case_experimental_ca_readback_launch_20260927.json');launch=json.loads(lp.read_text());sources[str(lp)]=sha(lp)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact complete CA readback',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0')
    assert sha(launch['cmdline'][1])==launch['script_sha256']
    audit_path=BASE/'experimental_structures/whole-domain-case-ca-readback-20260927-v1/receipt.json'
    audit=json.loads(audit_path.read_text());assert audit['status']=='complete_full_case_CA_raw_atom_and_grid_readback';sources[str(audit_path)]=sha(audit_path)
    mapping=BASE/'experimental_structures/whole-domain-case-ca-mapping-20260927-v1'
    cp=mapping/'config.json';assert audit['source_hashes'][str(cp)]==sha(cp);sources[str(cp)]=sha(cp);excluded=json.loads(cp.read_text())['deferred_entries']
    structural=BASE/'structural_comparisons';case=structural/'whole-domain-case-dossiers-20260927-v1'
    whole=rows(checked(case,'whole_reference_links.tsv'));domain=rows(checked(case,'domain_reference_links.tsv'))
    keyfields=['guide','family','gene_a','gene_b','reference_gene']
    whole_lookup={tuple(r[k] for k in keyfields):r['triad_id'] for r in whole}
    association=defaultdict(set)
    for r in domain:association[r['triad_key']].add(whole_lookup[tuple(r[k] for k in keyfields)])
    assert len(association)==26 and all(len(v)==1 for v in association.values())
    domains=defaultdict(list);needed=set()
    for line in checked(structural/'duplication-domain-common-residues-20260927-v1','triads.jsonl').open():
        r=json.loads(line)
        if r['triad_key'] in association:
            domains[next(iter(association[r['triad_key']]))].append(r)
            needed.update(r['interval_'+role] for role in ROLES)
    bounds={}
    for line in checked(structural/'duplication-domain-inputs-20260926-v1','inputs.jsonl').open():
        r=json.loads(line)
        if r['mask']=='full' and r['interval_id'] in needed:bounds[r['interval_id']]=r
    quartet_root=BASE/'experimental_structures/whole-domain-case-sequence-quartets-20260927-v1'
    dispositions=json.loads(checked(quartet_root,'triad_dispositions.json').read_text())
    with gzip.open(checked(quartet_root,'quartet_residue_maps.jsonl.gz'),'rt') as f:quartets=list(map(json.loads,f))
    models={(m,v) for r in quartets for m,v in zip(r['model_ids'],r['model_versions'])};inputs={}
    for name in ['duplication-alignment-inputs-20260926-v1','duplication-reference-alignment-inputs-20260926-v1']:
        for line in checked(structural/name,'inputs.jsonl').open():
            r=json.loads(line);key=(r['model_id'],r['version'])
            if key in models:
                ik=(*key,r['mask'])
                if ik in inputs:assert inputs[ik]['sequence']==r['sequence'] and inputs[ik]['original_positions']==r['original_positions']
                inputs[ik]=r
    coordinates={}
    for model in models:
        for mask in ['full','plddt70']:
            r=inputs[(*model,mask)]
            if r['status']=='ready':
                seq,xyz,confidence=load_pdb(r);sources[r['path']]=sha(r['path'])
                coordinates[(*model,mask)]={pos:xyz[i] for i,pos in enumerate(r['original_positions'])}
            else:coordinates[(*model,mask)]={}
    byentry=defaultdict(list)
    for r in quartets:byentry[r['entity_id'].rsplit('_',1)[0]].append(r)
    OUTPUT.mkdir(exist_ok=False);fit_path=OUTPUT/'pair_fits.tsv';map_path=OUTPUT/'observed_quartet_partitions.jsonl.gz'
    n_maps=n_fits=expected_maps=0;maxerror=0.;exclusions=[]
    with fit_path.open('w') as f,gzip.open(map_path,'wt') as mf:
        writer=None
        for entry,entry_quartets in sorted(byentry.items()):
            if entry in excluded:
                exclusions.append(dict(entry_id=entry,reason=excluded[entry],quartet_maps=len(entry_quartets)));continue
            path=mapping/(entry+'.residues.tsv.gz');assert sha(path)==audit['source_hashes'][str(path)];sources[str(path)]=sha(path)
            observed=defaultdict(dict)
            with gzip.open(path,'rt') as handle:
                for r in csv.DictReader(handle,delimiter='\t'):
                    key=(entry+'_'+r['entity_id'],r['model_number'],r['label_asym_id']);observed.setdefault(key,{})
                    if r['CA_status']=='unambiguous_full_occupancy_CA':
                        atom=json.loads(r['atom_records_json'])[0];observed[key][int(r['label_seq_id'])]=np.array([float(atom[k]) for k in ['Cartn_x','Cartn_y','Cartn_z']])
            for q in entry_quartets:
                modelkeys=list(zip(q['model_ids'],q['model_versions']))
                for index,key in enumerate(modelkeys):
                    full=inputs[(*key,'full')]
                    assert full['original_length']==q['query_lengths'][index]
                    assert 'S'+hashlib.sha256(full['sequence'].encode()).hexdigest()==q['sequence_ids'][index]
                chain_count=sum(entity==q['entity_id'] for entity,model,chain in observed)
                assert chain_count>0 and domains[q['triad_id']]
                expected_maps+=chain_count*2*len(domains[q['triad_id']])
                for (entity,exp_model,chain),ec in observed.items():
                    if entity!=q['entity_id']:continue
                    for mask in ['full','plddt70']:
                        coords=[coordinates[(*key,mask)] for key in modelkeys]+[ec]
                        common=[p for p in q['common_positions'] if all(p[i] in coords[i] for i in range(4))]
                        for dt in domains[q['triad_id']]:
                            bb=[bounds[dt['interval_'+role]] for role in ROLES]
                            assert [(b['model_id'],b['version']) for b in bb]==modelkeys
                            inside=[];outside=[];mixed=[]
                            for p in common:
                                flags=[b['start']<=p[i]<=b['end'] for i,b in enumerate(bb)]
                                (inside if all(flags) else outside if not any(flags) else mixed).append(p)
                            assert len(common)==len(inside)+len(outside)+len(mixed)
                            ident=dict(triad_id=q['triad_id'],domain_triad=dt['triad_key'],entity_id=entity,experimental_model=exp_model,label_asym_id=chain,mask=mask,context_indices_json=json.dumps(q['context_indices']))
                            record=dict(ident,common_positions=common,inside_positions=inside,outside_positions=outside,mixed_positions=mixed,original_query_lengths=q['query_lengths'],domain_lengths=[b['end']-b['start']+1 for b in bb],outside_lengths=[b['original_length']-(b['end']-b['start']+1) for b in bb])
                            mf.write(json.dumps(record,separators=(',',':'))+'\n');n_maps+=1
                            for region,positions in [('whole',common),('domain',inside),('outside',outside)]:
                                for pair,i,j in PAIRS:
                                    row=dict(ident,region=region,pair=pair,residues=len(positions),status='too_few_residues',rmsd='',relative_rotation_curvature='',outside_under_domain_transform_rmsd='')
                                    if len(positions)>=3:
                                        x,y=[np.array([coords[k][p[k]] for p in positions]) for k in [i,j]];g=geometry(x,y);row['relative_rotation_curvature']=g['relative_rotation_curvature']
                                        if g['geometry_status']!='unique_at_numeric_tolerance':row['status']='degenerate_rotation'
                                        else:
                                            rot,trans=transform(x,y,'svd');qr,qt=transform(x,y,'quaternion');value=residual(x,y,rot,trans);error=abs(value-residual(x,y,qr,qt));row.update(status='computed',rmsd=value)
                                            if region=='domain' and outside:
                                                ox,oy=[np.array([coords[k][p[k]] for p in outside]) for k in [i,j]]
                                                value=residual(ox,oy,rot,trans);row['outside_under_domain_transform_rmsd']=value;error=max(error,abs(value-residual(ox,oy,qr,qt)))
                                            assert error<1e-8;maxerror=max(maxerror,error)
                                    if writer is None:writer=csv.DictWriter(f,list(row),delimiter='\t',lineterminator='\n');writer.writeheader()
                                    writer.writerow(row);n_fits+=1
            print(entry,n_maps,n_fits,flush=True)
    assert n_maps==expected_maps and n_fits==18*n_maps
    with fit_path.open() as f:assert sum(1 for _ in csv.DictReader(f,delimiter='\t'))==n_fits
    with gzip.open(map_path,'rt') as f:assert sum(1 for _ in map(json.loads,f))==n_maps
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='complete_case_experimental_quartet_geometric_fits',source_hashes=sources,script_sha256=sha(__file__),fixtures=fixtures,partition_maps=n_maps,pair_fit_rows=n_fits,max_svd_quaternion_disagreement=maxerror,excluded_entries=exclusions,all_triad_dispositions=dispositions,artifacts={p.name:sha(p) for p in [fit_path,map_path]},scope='All common sequence hits, chain/models, masks and domain boundaries retained. Six pairwise proper-rotation fits use identical four-way observed residue sets per region. Domain transforms also evaluated outside without refitting. Degenerate and insufficient sets explicit. Coverage screens, full serialized numeric readback, construct dependence and scientific interpretation remain downstream. No independent validation or biological mechanism claim.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
