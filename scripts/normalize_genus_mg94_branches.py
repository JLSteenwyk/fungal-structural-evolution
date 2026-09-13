#!/usr/bin/env python3
"""Apply and independently verify upstream MG94 site-opportunity normalization."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
from itertools import product
import json
import math
from pathlib import Path
import re
import subprocess
import numpy as np
from Bio.Data import CodonTable
from audit_genus_codon_trees import sha


def opportunities(code):
    table=CodonTable.unambiguous_dna_by_id[code]
    codons=[''.join(c) for c in product('ACGT',repeat=3) if ''.join(c) not in table.stop_codons]
    counts=[]
    for codon in codons:
        s=n=0.
        for position in range(3):
            for base in 'ACGT':
                if base==codon[position]:continue
                mutant=codon[:position]+base+codon[position+1:]
                if mutant in table.stop_codons:continue
                if table.forward_table[mutant]==table.forward_table[codon]:s+=1/3
                else:n+=1/3
        counts.append((s,n))
    return codons,np.asarray(counts)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['fits','audit','binding','source','plan','install','opportunity-helper','output']:parser.add_argument('--'+name,type=Path,required=True)
    a=parser.parse_args()
    if a.output.exists():raise FileExistsError('Use a new output directory')
    plan=json.loads(a.plan.read_text());source=json.loads(a.source.read_text());fit=json.loads((a.fits/'receipt.json').read_text());audit=json.loads((a.audit/'receipt.json').read_text());binding=json.loads(a.binding.read_text())
    for path,key in [(a.fits/'receipt.json','source_fit_receipt_sha256'),(a.audit/'receipt.json','source_audit_sha256'),(a.source,'source_definition_sha256'),(a.binding,'binding_sha256')]:
        if sha(path)!=plan[key]:raise ValueError('Changed planned source')
    if audit['status']!='passed_saved_fit_readback' or audit['scope']!='all_executed_cases' or binding['source_fit_receipt_sha256']!=sha(a.fits/'receipt.json'):raise ValueError('Full audited/bound fits required')
    config=json.loads((a.fits/'config.json').read_text());exe=a.install/'bin/hyphy'
    if sha(exe)!=config['executable_sha256']:raise ValueError('Changed executable')
    helper=a.install/'share/hyphy/TemplateBatchFiles/libv3/tasks/genetic_code.bf'
    if sha(helper)!=source['installed_version_source_genetic_code_sha256']:raise ValueError('Changed original opportunity helper')
    if sha(a.opportunity_helper)!=plan['corrected_helper_sha256']:raise ValueError('Changed corrected opportunity helper')
    expected={r['case_id']:r['receipt_sha256'] for r in fit['case_receipts']}
    checked={r['case_id']:r['source_receipt_sha256'] for r in audit['source_cases']}
    if expected!=checked or len(expected)!=plan['cases']:raise ValueError('Source fit/audit grids differ')
    a.output.mkdir(parents=True)
    opp={code:opportunities(code) for code in [1,12]}
    prelude='skipCodeSelectionStep=1;\nLoadFunctionLibrary("TemplateModels/chooseGeneticCode.def");\nLoadFunctionLibrary('+json.dumps(str(a.opportunity_helper.resolve()))+');\n'
    test=prelude
    for code,index in [(1,0),(12,8)]:
        test+='ApplyGeneticCodeTable('+str(index)+');\ncounts=genetic_code.ComputePairwiseDifferencesAndExpectedSites(_Genetic_Code,{});\n'
        for i,(codon,values) in enumerate(zip(*opp[code])):
            for name,value in zip(['SS','NS'],values):test+='assert(Abs((counts[genetic_code.'+name+'])['+str(i)+']-('+repr(float(value))+'))<1e-12,"Opportunity mismatch code '+str(code)+' '+codon+' '+name+'");\n'
    test+='fprintf(stdout,"ALL_OPPORTUNITIES_PASSED\\n");\n'
    (a.output/'check_opportunities.bf').write_text(test)
    with (a.output/'check_opportunities.log').open('w') as stream:subprocess.run([str(exe.resolve()),'CPU=1',str((a.output/'check_opportunities.bf').resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True)
    def run(case):
        folder=a.fits/case;rp=folder/'receipt.json'
        if sha(rp)!=expected[case]:raise ValueError('Changed completed fit')
        r=json.loads(rp.read_text());rc=json.loads((folder/'config.json').read_text())
        for name,digest in r['artifacts'].items():
            if sha(folder/name)!=digest:raise ValueError('Changed fitted artifact')
        code=int(rc['translation_table']);codons,counts=opp[code]
        saved=(folder/'fit.bf').read_text()
        match=re.findall(r'^([A-Za-z0-9_.]+\.model_MGREV_pi)=\{\s*\n(.*?)\n\};',saved,re.M|re.S)
        if len(match)!=1:raise ValueError('Ambiguous equilibrium frequency vector')
        variable,body=match[0];tokens=re.findall(r'\{([^{}]+)\}',body)
        freq=np.array([float(x.strip()) for x in tokens])
        if len(freq)!=len(codons) or not np.isfinite(freq).all() or (freq<0).any() or abs(freq.sum()-1)>1e-10:raise ValueError('Invalid frequency vector')
        s,n=freq@counts
        if min(s,n)<=0:raise ValueError('Nonpositive weighted opportunity')
        output=a.output/case;output.mkdir()
        hbl=prelude+'ApplyGeneticCodeTable('+str(0 if code==1 else 8)+');\ncounts=genetic_code.ComputePairwiseDifferencesAndExpectedSites(_Genetic_Code,{});\nExecuteAFile('+json.dumps(str((folder/'fit.bf').resolve()))+');\nS=+('+variable+' $ counts[genetic_code.SS]);\nNS=+('+variable+' $ counts[genetic_code.NS]);\nfprintf(stdout,"EXPECTED_S=",Format(S,0,16),"\\nEXPECTED_NS=",Format(NS,0,16),"\\nFREQ_SUM=",Format(+'+variable+',0,16),"\\n");\n'
        (output/'check.bf').write_text(hbl)
        with (output/'check.log').open('w') as stream:subprocess.run([str(exe.resolve()),'CPU=1',str((output/'check.bf').resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True)
        log=(output/'check.log').read_text();values={k:float(v) for k,v in re.findall(r'^(EXPECTED_S|EXPECTED_NS|FREQ_SUM)=(\S+)',log,re.M)}
        if set(values)!={'EXPECTED_S','EXPECTED_NS','FREQ_SUM'} or not np.allclose([s,n,freq.sum()],[values['EXPECTED_S'],values['EXPECTED_NS'],values['FREQ_SUM']],rtol=0,atol=1e-12):raise ValueError('Independent weighted opportunities differ')
        report=json.loads((folder/'fit.json').read_text());branchrows=[]
        for node,row in report['branch attributes']['0'].items():
            ds=row['synonymous']*(s+n)/s;dn=row['nonsynonymous']*(s+n)/n
            if not math.isfinite(ds+dn) or min(ds,dn)<0:raise ValueError('Invalid normalized distance')
            branchrows.append({'case_id':case,'node':node,'translation_table':code,'synonymous_component_per_site':row['synonymous'],'nonsynonymous_component_per_site':row['nonsynonymous'],'dS_upstream_equal_alternative_convention':ds,'dN_upstream_equal_alternative_convention':dn,'normalized_dN_over_dS':dn/ds if ds>0 else '', 'fitted_global_omega':r['reported_omega']})
        summary={'case_id':case,'translation_table':code,'sense_codons':len(codons),'frequency_sum':float(freq.sum()),'expected_synonymous_sites':float(s),'expected_nonsynonymous_sites':float(n),'synonymous_scale':float((s+n)/s),'nonsynonymous_scale':float((s+n)/n),'branches':len(branchrows),'maximum_branch_dS':max(x['dS_upstream_equal_alternative_convention'] for x in branchrows),'maximum_branch_dN':max(x['dN_upstream_equal_alternative_convention'] for x in branchrows),'independent_opportunity_max_difference':max(abs(s-values['EXPECTED_S']),abs(n-values['EXPECTED_NS'])),'marker_copy_caveat':rc['marker_copy_caveat']}
        proof={'case_id':case,'source_fit_receipt_sha256':sha(rp),'check_script_sha256':sha(output/'check.bf'),'check_log_sha256':sha(output/'check.log')}
        return summary,branchrows,proof
    summaries=[];branches=[];proofs=[]
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        futures=[pool.submit(run,case) for case in sorted(expected)]
        for i,future in enumerate(as_completed(futures),1):
            summary,rows,proof=future.result();summaries.append(summary);branches.extend(rows);proofs.append(proof)
            if i%100==0:print('Normalized',i,flush=True)
    if len(branches)!=plan['branches']:raise ValueError('Branch grid differs')
    summaries.sort(key=lambda x:x['case_id']);branches.sort(key=lambda x:(x['case_id'],x['node']))
    for name,rows in [('case_normalization.tsv',summaries),('normalized_branches.tsv',branches)]:
        with (a.output/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    result={'status':'complete_independently_checked_mg94_opportunity_normalization','cases':len(summaries),'branches':len(branches),'sense_codon_opportunity_checks':244,'per_fit_independent_checks':len(proofs),'source_fit_receipt_sha256':sha(a.fits/'receipt.json'),'source_audit_sha256':sha(a.audit/'receipt.json'),'source_definition_sha256':sha(a.source),'plan_sha256':sha(a.plan),'script_sha256':sha(Path(__file__)),'checks':sorted(proofs,key=lambda x:x['case_id']),'artifacts':{name:sha(a.output/name) for name in ['check_opportunities.bf','check_opportunities.log','case_normalization.tsv','normalized_branches.tsv']},'interpretation':plan['interpretation']}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['checks','artifacts']},indent=2))


if __name__=='__main__':main()
