#!/usr/bin/env python3
"""Check generated model semantics and finite initial probabilities, not convergence."""
import argparse, collections, json, math, re
from pathlib import Path
from Bio import Phylo
from prepare_case_ancestral_neighborhoods import sha


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',required=True,type=Path);a=p.parse_args()
    pp=Path('metadata/baliphy_initialization_plan_20260927.json');plan=json.loads(pp.read_text())
    for path,h in plan['pins'].items():assert sha(path)==h,path
    source=json.loads(Path(plan['input_receipt']).read_text());root=Path(plan['output'])
    out=a.output;out.mkdir(parents=True,exist_ok=False);rows=[];pins={str(pp):sha(pp),str(Path(__file__)):sha(__file__)}
    for job in source['jobs']:
        folder=root/job['job_id'];rp=folder/'receipt.json';row=dict(job_id=job['job_id'],status='pending_receipt')
        if not rp.exists():rows.append(row);continue
        r=json.loads(rp.read_text());assert r['plan_sha256']==sha(pp) and r['job']==job
        pins[str(rp)]=sha(rp)
        for f,h in r['artifacts'].items():assert sha(folder/f)==h,(folder,f)
        row.update(exit_code=r['exit_code'],producer_status=r['status'])
        if r['exit_code']!=0:
            row['status']='unsuccessful_initialization_retained';rows.append(row);continue
        code=(folder/'BAliPhy.Main.hs').read_text()
        for fragment in ['sample (symmetricDirichletOn (letterSet alpha) 1)','IModel.rs07 0.01 3 topology','SModel.gammaRatesOn 1 4','sample_scale  = 1','mkUnalignedCharacterData aa','dropInternalLabels','if isTest then printInitialModel']:
            assert fragment in code,(job['job_id'],fragment)
        assert str(Path(job['alignment']).resolve()) in code and str(Path(job['tree']).resolve()) in code
        text=re.sub(r'\x1b\[[0-9;]*m','',(folder/'stdout.log').read_text())
        line=[x for x in text.splitlines() if re.search(r'\biter = 0\b',x)];assert len(line)==1
        def number(name):
            match=re.search(r'(?:^|\s)'+re.escape(name)+r' = ([-+0-9.eE]+)(?:\s|$)',line[0]);assert match,(job['job_id'],name)
            value=float(match.group(1));assert math.isfinite(value);return value
        prior,likelihood,posterior=[number(n) for n in ['prior','likelihood','posterior']]
        assert abs(prior+likelihood-posterior)<1e-7
        freq={aa:float(value) for aa,value in re.findall(r'F:pi\[([A-Z])\] = ([-+0-9.eE]+)',line[0])}
        assert set(freq)==set('ARNDCQEGHILKMFPSTWYV') and all(0<v<1 for v in freq.values()) and abs(sum(freq.values())-1)<1e-10
        tree=Phylo.read(job['tree'],'newick');total=sum(n.branch_length or 0 for n in tree.find_clades() if n is not tree.root)
        assert abs(number('|T|')-total)<1e-9 and number('scale1')==1
        assert number('P1/|A|')>0
        row.update(status='initial_model_semantics_and_finite_scores_checked',reported_prior=prior,reported_likelihood=likelihood,reported_posterior=posterior,initial_alignment_columns=number('P1/|A|'),tree_length=number('|T|'),frequencies=freq,elapsed_seconds=r['elapsed_seconds'],peak_sampled_rss_bytes=r['peak_sampled_rss_bytes'])
        rows.append(row)
    (out/'dispositions.json').write_text(json.dumps(rows,indent=2)+'\n')
    counts=dict(collections.Counter(r['status'] for r in rows))
    receipt=dict(status='partial_initialization_snapshot' if counts.get('pending_receipt') else 'all_initialization_dispositions_audited',expected_jobs=324,dispositions=counts,pins=pins,artifacts={'dispositions.json':sha(out/'dispositions.json')},scope='Generated model code, fixed input paths, initial score arithmetic, frequency normalization and total tree length checked. No independent likelihood replay, runtime tip/node reconstruction readback, MCMC mixing, effective sample size or posterior qualification. All324 outcomes retained.')
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(counts))
if __name__=='__main__':main()
