#!/usr/bin/env python3
"""Export all coverage screens and paired fungal/experimental reference contrasts."""
import csv,gzip,json,subprocess,time
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import psutil
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results/experimental_structures')
ROOT=BASE/'whole-domain-case-quartet-fits-20260927-v1'
AUDIT=BASE/'whole-domain-case-quartet-readback-20260927-v1'
OUTPUT=BASE/'whole-domain-case-quartet-contrasts-20260927-v1'
ID=['triad_id','domain_triad','entity_id','experimental_model','label_asym_id','mask','context_indices_json']


def qualifies(count,lengths,n,c):
    direct=count>=n and all(length>0 and count*100>=c*length for length in lengths)
    rational=count>=n and all(length>0 and Fraction(count,length)>=Fraction(c,100) for length in lengths)
    assert direct==rational
    return int(direct)


def difference(left,right):
    value=float(left)-float(right)
    assert abs(Decimal(str(value))-(Decimal(left)-Decimal(right)))<Decimal('1e-10')
    return value


def main():
    lp=Path('metadata/case_experimental_quartet_readback_launch_20260927.json');launch=json.loads(lp.read_text());lh=sha(lp)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact full quartet readback',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0') and sha(launch['cmdline'][1])==launch['script_sha256']
    audit=json.loads((AUDIT/'receipt.json').read_text());assert audit['status']=='complete_full_experimental_quartet_readback'
    r=json.loads((ROOT/'receipt.json').read_text());assert audit['source_hashes'][str(ROOT/'receipt.json')]==sha(ROOT/'receipt.json')
    sources={str(lp):lh,str(AUDIT/'receipt.json'):sha(AUDIT/'receipt.json'),str(ROOT/'receipt.json'):sha(ROOT/'receipt.json')}
    for name,digest in r['artifacts'].items():assert sha(ROOT/name)==digest;sources[str(ROOT/name)]=digest
    OUTPUT.mkdir(exist_ok=False);sp=OUTPUT/'coverage_screens.tsv';cp=OUTPUT/'paired_reference_contrasts.tsv';screens=contrasts=partitions=0
    with (ROOT/'pair_fits.tsv').open() as f,gzip.open(ROOT/'observed_quartet_partitions.jsonl.gz','rt') as mf,sp.open('w') as sf,cp.open('w') as cf:
        fits=iter(csv.DictReader(f,delimiter='\t'));sw=cw=None
        for m in map(json.loads,mf):
            ident={k:m[k] for k in ID};partitions+=1;current={}
            for _ in range(18):
                row=next(fits);assert all(row[k]==m[k] for k in ID);key=(row['region'],row['pair']);assert key not in current;current[key]=row
            assert set(current)=={(region,pair) for region in ['whole','domain','outside'] for pair in ['ab','ar','br','ae','be','re']}
            counts=[len(m[k]) for k in ['common_positions','inside_positions','outside_positions']]
            lengths=[m[k] for k in ['original_query_lengths','domain_lengths','outside_lengths']]
            for n in [30,50]:
                for c in [50,70,90]:
                    whole,domain,outside=[qualifies(count,denom,n,c) for count,denom in zip(counts,lengths)]
                    row=dict(ident,screen=f'n{n}_c{c}',whole_residues=counts[0],domain_residues=counts[1],outside_residues=counts[2],mixed_residues=len(m['mixed_positions']),whole_pass=whole,domain_pass=domain,outside_pass=outside,domain_and_outside_pass=domain*outside)
                    if sw is None:sw=csv.DictWriter(sf,list(row),delimiter='\t',lineterminator='\n');sw.writeheader()
                    sw.writerow(row);screens+=1
            for label,region,metric in [('whole','whole','rmsd'),('domain','domain','rmsd'),('outside_independent','outside','rmsd'),('outside_domain_anchored','domain','outside_under_domain_transform_rmsd')]:
                row=dict(ident,metric=label,coverage_requirement='domain_and_outside' if label=='outside_domain_anchored' else 'outside' if label=='outside_independent' else label,
                         fungal_reference_status='unavailable_fit',experimental_reference_status='unavailable_fit',fungal_reference_ar_minus_br='',experimental_reference_ae_minus_be='',experimental_minus_fungal_contrast='')
                for prefix,left,right,column in [('fungal','ar','br','fungal_reference_ar_minus_br'),('experimental','ae','be','experimental_reference_ae_minus_be')]:
                    a,b=current[region,left],current[region,right]
                    if a['status']=='computed' and b['status']=='computed' and a[metric]!='' and b[metric]!='':
                        row[prefix+'_reference_status']='computed';row[column]=difference(a[metric],b[metric])
                    row[prefix+'_left_fit_status']=a['status'];row[prefix+'_right_fit_status']=b['status']
                if row['fungal_reference_status']==row['experimental_reference_status']=='computed':
                    row['experimental_minus_fungal_contrast']=difference(str(row['experimental_reference_ae_minus_be']),str(row['fungal_reference_ar_minus_br']))
                if cw is None:cw=csv.DictWriter(cf,list(row),delimiter='\t',lineterminator='\n');cw.writeheader()
                cw.writerow(row);contrasts+=1
        assert next(fits,None) is None
    assert partitions==r['partition_maps'] and screens==6*partitions and contrasts==4*partitions
    for path,count in [(sp,screens),(cp,contrasts)]:
        with path.open() as f:assert sum(1 for _ in csv.DictReader(f,delimiter='\t'))==count
    for path,digest in sources.items():assert sha(path)==digest
    result=dict(status='complete_case_quartet_coverage_and_paired_contrasts',source_hashes=sources,script_sha256=sha(__file__),partition_maps=partitions,coverage_rows=screens,contrast_rows=contrasts,all_triad_dispositions=r['all_triad_dispositions'],artifacts={p.name:sha(p) for p in [sp,cp]},scope='All six original-length coverage screens and four paired contrasts per partition retained, including unavailable fits and below-threshold cases. AR-BR uses fungal reference; AE-BE uses the experimental homolog on the same residues. Differences are descriptive and reference-dependent, not an evolutionary polarity, independent replication, prediction-accuracy estimate, effect uncertainty or selection test. Cases were selected for prior whole/domain contrast; no confirmatory p-values.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
