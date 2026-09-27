#!/usr/bin/env python3
"""Independently check every quartet coverage row and signed contrast export."""
import csv,gzip,json,subprocess,time
from collections import Counter
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
import psutil
from screen_duplication_domain_alignment_coverage import sha

BASE=Path('results/experimental_structures')
ROOT=BASE/'whole-domain-case-quartet-contrasts-20260927-v1'
FITS=BASE/'whole-domain-case-quartet-fits-20260927-v1'
OUTPUT=BASE/'whole-domain-case-quartet-contrast-readback-20260927-v1'
ID=['triad_id','domain_triad','entity_id','experimental_model','label_asym_id','mask','context_indices_json']
METRICS={'whole':('whole','rmsd','whole'),'domain':('domain','rmsd','domain'),'outside_independent':('outside','rmsd','outside'),'outside_domain_anchored':('domain','outside_under_domain_transform_rmsd','domain_and_outside')}


def rows(p):
    with p.open() as f:yield from csv.DictReader(f,delimiter='\t')


def main():
    lp=Path('metadata/case_experimental_quartet_contrasts_launch_20260927.json');launch=json.loads(lp.read_text());lh=sha(lp)
    while True:
        try:
            p=psutil.Process(launch['pid'])
            if p.create_time()!=launch['created'] or p.status()==psutil.STATUS_ZOMBIE:break
            assert p.cmdline()==launch['cmdline']
        except psutil.NoSuchProcess:break
        print('Waiting for exact contrast exporter',launch['pid'],flush=True);time.sleep(30)
    state=dict(x.split('=',1) for x in subprocess.check_output(['systemctl','--user','show',launch['unit'],'-p','ActiveState','-p','Result','-p','ExecMainStatus'],text=True).splitlines())
    assert state==dict(ActiveState='inactive',Result='success',ExecMainStatus='0') and sha(launch['cmdline'][1])==launch['script_sha256']
    rp=ROOT/'receipt.json';r=json.loads(rp.read_text());assert r['status']=='complete_case_quartet_coverage_and_paired_contrasts'
    sources={**r['source_hashes'],str(rp):sha(rp),str(lp):lh}
    for name,digest in r['artifacts'].items():sources[str(ROOT/name)]=digest
    for p,digest in sources.items():assert sha(p)==digest
    partitions={}
    with gzip.open(FITS/'observed_quartet_partitions.jsonl.gz','rt') as f:
        for m in map(json.loads,f):
            key=tuple(m[k] for k in ID);assert key not in partitions;partitions[key]=m
    assert len(partitions)==r['partition_maps']
    fits={}
    for row in rows(FITS/'pair_fits.tsv'):
        key=tuple(row[k] for k in ID)+(row['region'],row['pair']);assert key not in fits;fits[key]=row
    assert len(fits)==18*len(partitions)
    expected={(key,f'n{n}_c{c}') for key in partitions for n in [30,50] for c in [50,70,90]};screen_count=0
    for row in rows(ROOT/'coverage_screens.tsv'):
        key=tuple(row[k] for k in ID);ek=(key,row['screen']);assert ek in expected;expected.remove(ek);screen_count+=1;m=partitions[key]
        n,c=(int(x[1:]) for x in row['screen'].split('_'));flags={}
        for name,position_field,length_field in [('whole','common_positions','original_query_lengths'),('domain','inside_positions','domain_lengths'),('outside','outside_positions','outside_lengths')]:
            count=len(m[position_field]);assert int(row[name+'_residues'])==count
            flag=count>=n and all(length>0 and Fraction(count,length)>=Fraction(c,100) for length in m[length_field]);assert int(row[name+'_pass'])==flag;flags[name]=flag
        assert int(row['mixed_residues'])==len(m['mixed_positions']) and int(row['domain_and_outside_pass'])==(flags['domain'] and flags['outside'])
    assert not expected and screen_count==r['coverage_rows']
    expected={(key,metric) for key in partitions for metric in METRICS};contrast_count=0;counts=Counter();maximum=0.
    for row in rows(ROOT/'paired_reference_contrasts.tsv'):
        key=tuple(row[k] for k in ID);ek=(key,row['metric']);assert ek in expected;expected.remove(ek);contrast_count+=1
        region,field,requirement=METRICS[row['metric']];assert row['coverage_requirement']==requirement;values={}
        for prefix,a,b,column in [('fungal','ar','br','fungal_reference_ar_minus_br'),('experimental','ae','be','experimental_reference_ae_minus_be')]:
            left,right=fits[(*key,region,a)],fits[(*key,region,b)]
            assert row[prefix+'_left_fit_status']==left['status'] and row[prefix+'_right_fit_status']==right['status']
            available=left['status']==right['status']=='computed' and left[field]!='' and right[field]!=''
            assert row[prefix+'_reference_status']==('computed' if available else 'unavailable_fit')
            if available:
                value=Decimal(left[field])-Decimal(right[field]);error=abs(Decimal(row[column])-value);assert error<Decimal('1e-10');maximum=max(maximum,float(error));values[prefix]=value
            else:assert row[column]==''
        if len(values)==2:
            error=abs(Decimal(row['experimental_minus_fungal_contrast'])-(values['experimental']-values['fungal']));assert error<Decimal('1e-10');maximum=max(maximum,float(error));counts[row['metric']+':both_available']+=1
        else:assert row['experimental_minus_fungal_contrast']=='';counts[row['metric']+':unavailable']+=1
    assert not expected and contrast_count==r['contrast_rows']
    for p,digest in sources.items():assert sha(p)==digest
    OUTPUT.mkdir(exist_ok=False);result=dict(status='complete_full_case_quartet_contrast_readback',terminal_state=state,source_hashes=sources,script_sha256=sha(__file__),partition_maps=len(partitions),coverage_rows=screen_count,contrast_rows=contrast_count,availability_counts=dict(counts),maximum_decimal_subtraction_error=maximum,scope='All six screens and all four metrics verified for every partition. Counts reconstructed from source positions, rational coverage thresholds and original denominators; signed contrasts and all unavailable statuses rechecked against source pair fits using Decimal arithmetic. Descriptive selected-case outputs only; no inference of biological effect or independence.')
    (OUTPUT/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
