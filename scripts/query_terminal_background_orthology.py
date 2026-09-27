#!/usr/bin/env python3
"""Query both audited native ortholog streams for every modeled background candidate."""
import argparse,csv,json,subprocess,time
from pathlib import Path
import psutil
from run_ortholog_pair_guide_comparison import sha


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--plan',type=Path,required=True);a=ap.parse_args()
    plan=json.loads(a.plan.read_text());ph=sha(a.plan)
    def verify():
        assert sha(a.plan)==ph
        for path,h in plan['pins'].items():assert sha(path)==h,path
    verify();dep=plan['dependency']
    while True:
        try:
            p=psutil.Process(dep['pid'])
            if p.create_time()!=dep['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==dep['cmdline']
        except psutil.NoSuchProcess:break
        time.sleep(30)
    verify();source=Path(plan['source']);receipt=json.loads((source/'receipt.json').read_text());rh=sha(source/'receipt.json')
    assert receipt['status']=='complete_terminal_sister_guide_comparison_pending_readback'
    table=source/'modeled_candidate_union.tsv';assert sha(table)==receipt['artifacts'][table.name]
    with table.open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==receipt['modeled_candidate_union']
    genes={r[k] for r in rows for k in ['gene_a','gene_b']}
    ids={};sp={}
    for line in Path(plan['species_ids']).read_text().splitlines():
        native,label=line.split(': ',1);sp[native]=label.rsplit('.',1)[0]
    with Path(plan['sequence_ids']).open() as f:
        for ordinal,line in enumerate(f):
            native,protein=line.rstrip('\n').split(': ',1);gene=sp[native.split('_')[0]]+'_'+protein
            if gene in genes:
                assert gene not in ids and ordinal<2**24
                ids[gene]=ordinal
    assert set(ids)==genes
    keyed={}
    for row in rows:
        x,y=sorted([ids[row['gene_a']],ids[row['gene_b']]]);assert x<y
        key=f'{x:06x}{y:06x}';assert key not in keyed;keyed[key]=row
    out=Path(plan['output']);out.mkdir(exist_ok=False)
    queries=out/'queries.hex';queries.write_text(''.join(k+'\n' for k in sorted(keyed)))
    binary=out/'query_membership'
    subprocess.run([plan['compiler'],'-O2','-std=c++17',plan['cpp'],'-o',str(binary)],check=True)
    subprocess.run([plan['python'],plan['fixtures'],'--binary',str(binary)],check=True)
    stats={};members={}
    for source_stream in plan['streams']:
        guide=source_stream['guide'];r=json.loads(Path(source_stream['receipt']).read_text())
        assert r['status']=='complete_native_ortholog_pair_multiplicity_audit'
        for field in ['duplicate_directed_incidences','pairs_missing_reverse','pairs_with_repeated_direction','pairs_with_unequal_multiplicity']:assert r[field]==0
        assert sha(source_stream['stream'])==r['sorted_pairs_sha256']
        result=out/(guide+'_membership.tsv')
        run=subprocess.run([str(binary),source_stream['stream'],str(queries),str(result)],capture_output=True,text=True,check=True)
        counts=json.loads(run.stdout);assert counts['queries']==len(rows) and counts['stream_pairs']==r['unique_unordered_pairs']
        stats[guide]=counts;members[guide]={}
        with result.open() as f:
            for line in f:
                key,present=line.rstrip().split('\t');assert key not in members[guide] and present in ['0','1'];members[guide][key]=present
        assert set(members[guide])==set(keyed)
        assert sum(int(v) for v in members[guide].values())==counts['present']
        print(guide,json.dumps(counts),flush=True)
    with (out/'candidate_orthology_membership.tsv').open('w') as f:
        fields=list(rows[0])+['native_pair_key','profile_native_ortholog','mafft_native_ortholog']
        writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader()
        for key,row in sorted(keyed.items()):writer.writerow(dict(row,native_pair_key=key,profile_native_ortholog=members['profile'][key],mafft_native_ortholog=members['mafft'][key]))
    verify();assert sha(source/'receipt.json')==rh and sha(table)==receipt['artifacts'][table.name]
    for s in plan['streams']:assert sha(s['stream'])==json.loads(Path(s['receipt']).read_text())['sorted_pairs_sha256']
    result=dict(status='complete_background_native_ortholog_membership_pending_readback',plan_sha256=ph,source_receipt_sha256=rh,candidate_rows=len(rows),guides=stats,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Exact native output membership for every modeled candidate in either guide, using source-line protein ordinals and full reciprocal-stream scans. All absence and disagreement retained. Source guide comparison and this membership stage still require independent readback. Native orthology assignment is not independent biological validation or matched-control eligibility.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
