#!/usr/bin/env python3
"""Map all configurations to exact effective-input groups without deleting labels."""
import collections,hashlib,json,re
from pathlib import Path
from Bio import SeqIO
from prepare_case_ancestral_neighborhoods import sha


def digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True,separators=(',',':')).encode()).hexdigest()


def main():
    pp=Path('metadata/baliphy_initialization_plan_20260927.json');plan=json.loads(pp.read_text())
    for p,h in plan['pins'].items():assert sha(p)==h,p
    ip=Path(plan['input_receipt']);source=json.loads(ip.read_text());groups=collections.defaultdict(list);rows=[];pins={str(pp):sha(pp),str(ip):sha(ip),str(Path(__file__)):sha(__file__)}
    library=Path('data/software_audits/baliphy-4.3-20260927/install/bali-phy-4.3/lib/bali-phy/haskell/Bio/Sequence.hs')
    assert "indices' = map (\\(label,is) -> (label, stripGaps is)) indices" in library.read_text()
    pins[str(library)]=sha(library)
    evidence=collections.defaultdict(list)
    for j in source['jobs']:
        assert sha(j['alignment'])==j['alignment_sha256'] and sha(j['tree'])==j['tree_sha256']
        rs=list(SeqIO.parse(j['alignment'],'fasta'))
        assert len(rs)==len({r.id for r in rs})==j['proteins']
        assert all(r.description==r.id and set(str(r.seq))<=set('ARNDCQEGHILKMFPSTWYVX-') for r in rs)
        ordered=[(r.id,str(r.seq).replace('-','')) for r in rs]
        payload=dict(ordered_ungapped_sequences=ordered,tree_bytes_sha256=sha(j['tree']),binary_sha256=sha(plan['binary']),options=plan['options'])
        key=digest(payload);groups[key].append(j['job_id'])
        rows.append(dict(job_id=j['job_id'],input_id=j['input_id'],group=key,proteins=j['proteins'],alignment_sha256=j['alignment_sha256'],tree_sha256=j['tree_sha256'],ordered_sequence_sha256=digest(ordered)))
        rp=Path(plan['output'])/j['job_id']/'receipt.json'
        if not rp.exists():continue
        r=json.loads(rp.read_text());assert r['job']==j and r['plan_sha256']==sha(pp);pins[str(rp)]=sha(rp)
        for f,h in r['artifacts'].items():assert sha(rp.parent/f)==h
        if r['exit_code']!=0:continue
        code=(rp.parent/'BAliPhy.Main.hs').read_text()
        for path,token in [(j['alignment'],'__ALIGNMENT__'),(j['tree'],'__TREE__')]:
            literal=json.dumps(str(Path(path).resolve()));assert code.count(literal)==1;code=code.replace(literal,json.dumps(token))
        text=re.sub(r'\x1b\[[0-9;]*m','',(rp.parent/'stdout.log').read_text())
        lines=[re.sub(r'\s+',' ',line).strip() for line in text.splitlines() if re.search(r'\biter = 0\b',line)]
        assert len(lines)==1
        evidence[key].append(dict(job_id=j['job_id'],normalized_program_sha256=hashlib.sha256(code.encode()).hexdigest(),initial_state_sha256=hashlib.sha256(lines[0].encode()).hexdigest()))
    assert len(rows)==324 and len({r['job_id'] for r in rows})==324
    comparisons=[]
    for key,values in evidence.items():
        if len(values)>1:
            comparisons.append(dict(group=key,completed_members=len(values),normalized_program_variants=len({v['normalized_program_sha256'] for v in values}),initial_state_variants=len({v['initial_state_sha256'] for v in values})))
    out=Path('results/ancestral/baliphy-input-equivalence-20260927-v1');out.mkdir(parents=True,exist_ok=False)
    (out/'configuration_mapping.json').write_text(json.dumps(rows,indent=2)+'\n')
    (out/'completed_equivalence_checks.json').write_text(json.dumps(dict(evidence=evidence,comparisons=comparisons),indent=2)+'\n')
    result=dict(status='all_configurations_mapped_observed_initialization_equivalence_assessed',configurations=len(rows),exact_effective_input_groups=len(groups),group_size_counts=dict(collections.Counter(len(v) for v in groups.values())),completed_groups_compared=len(comparisons),groups_with_program_differences=sum(v['normalized_program_variants']>1 for v in comparisons),groups_with_initial_state_differences=sum(v['initial_state_variants']>1 for v in comparisons),pins=pins,artifacts={p.name:sha(p) for p in out.iterdir()},scope='Exact ordered sequences, complete tree bytes, binary and model/seed options; strip only gaps as pinned loader does. All proteins, roots, resolutions, floors and original labels retained. No existing job cancelled or replaced. Does not merge different trees, seeds, priors, models or genuine independent chains; future reuse must retain all label mappings and compare model code. Initial-state agreement is not MCMC convergence or posterior equivalence validation.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['pins','artifacts','scope']}))
if __name__=='__main__':main()
