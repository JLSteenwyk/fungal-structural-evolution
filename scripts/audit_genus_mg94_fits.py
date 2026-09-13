#!/usr/bin/env python3
"""Audit MG94 outputs and reevaluate every saved likelihood without optimization."""
import argparse
import csv
import io
import json
import math
from pathlib import Path
import re
import subprocess
from Bio import Phylo, SeqIO
from audit_genus_codon_trees import sha, split_map


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for key in ['fits','install','output']:
        parser.add_argument('--'+key,required=True,type=Path)
    parser.add_argument('--allow-incomplete',action='store_true')
    args=parser.parse_args()
    if args.output.exists():raise FileExistsError('Use a new immutable audit directory')
    cp=args.fits/'config.json';config=json.loads(cp.read_text())
    paths=sorted(args.fits.glob('*/receipt.json'))
    if not paths:raise ValueError('No completed fits')
    if not args.allow_incomplete and len(paths)!=config['cases']:raise ValueError('Full execution is incomplete')
    if len(paths)>config['cases']:raise ValueError('Too many cases')
    exe=args.install/'bin/hyphy'
    if sha(exe)!=config['executable_sha256']:raise ValueError('Changed executable')
    args.output.mkdir(parents=True)
    summaries=[];branches=[];sources=[]
    for count,rp in enumerate(paths,1):
        folder=rp.parent;receipt=json.loads(rp.read_text());case=folder.name
        rc=json.loads((folder/'config.json').read_text())
        if receipt['status']!='complete_mg94_execution_pending_independent_audit' or receipt['case_id']!=case or receipt['config_sha256']!=sha(folder/'config.json') or rc['parent_config_sha256']!=sha(cp):raise ValueError('Changed case provenance')
        for name,digest in receipt['artifacts'].items():
            if sha(folder/name)!=digest:raise ValueError('Changed fit artifact')
        command=rc['command']
        for option,value in {'--type':'global','--frequencies':'CF3x4','--lrt':'No','--kill-zero-lengths':'No','--code':{'1':'Universal','12':'Alt-Yeast-Nuclear'}[rc['translation_table']]}.items():
            if command[command.index(option)+1]!=value:raise ValueError('Wrong model or code setting')
        alignment=Path(command[command.index('--alignment')+1])
        if sha(alignment)!=rc['alignment_sha256'] or sha(folder/'tree.nwk')!=rc['tree_sha256']:raise ValueError('Changed model inputs')
        records=list(SeqIO.parse(alignment,'fasta'));taxa={r.id for r in records};ncols=len(records[0].seq)//3
        result=json.loads((folder/'fit.json').read_text());fit=result['fits']['Standard MG94']
        if len(taxa)!=len(records) or result['input']['number of sequences']!=len(taxa) or result['input']['number of sites']!=ncols or result['input']['partition count']!=1:raise ValueError('Fit dimensions differ')
        if result['data partitions']['0']['coverage']!=[list(range(ncols))]:raise ValueError('Codon column coverage differs')
        original=Phylo.read(folder/'tree.nwk','newick');exported=Phylo.read(io.StringIO(result['input']['trees']['0']+';'),'newick')
        if set(split_map(original,taxa))!=set(split_map(exported,taxa)):raise ValueError('Exported topology differs')
        nodes={n.name for n in original.find_clades() if n is not original.root}
        attr=result['branch attributes']['0']
        if set(attr)!=nodes or len(attr)!=2*len(taxa)-3:raise ValueError('Branch identity grid differs')
        largest_component_error=0
        for node,values in attr.items():
            vals=[values[k] for k in ['Standard MG94','synonymous','nonsynonymous']]
            if any(not math.isfinite(v) or v<0 for v in vals):raise ValueError('Invalid branch component')
            error=abs(vals[0]-vals[1]-vals[2]);largest_component_error=max(error,largest_component_error)
            branches.append({'case_id':case,'node':node,'total_substitutions_per_site':vals[0],'synonymous_component_per_site':vals[1],'nonsynonymous_component_per_site':vals[2],'absolute_additivity_error':error})
        omega=fit['Rate Distributions']['non-synonymous/synonymous rate ratio'];ci=fit['Confidence Intervals']['non-synonymous/synonymous rate ratio'];ll=fit['Log Likelihood']
        if any(not math.isfinite(v) for v in [omega,ci['LB'],ci['UB'],ll]) or omega<0 or ci['LB']<0 or ci['UB']<ci['LB']:raise ValueError('Invalid estimate or interval')
        saved=(folder/'fit.bf').read_text()
        lfs=re.findall(r'^LikelihoodFunction\s+(\S+)\s*=',saved,re.M)
        omegas=re.findall(r'^global\s+[^\s=]+\.omega=([^;]+);',saved,re.M)
        if len(lfs)!=1 or len(omegas)!=1:raise ValueError('Saved model identity is ambiguous')
        saved_omega=float(omegas[0])
        check=args.output/case;check.mkdir()
        bf=check/'reevaluate.bf'
        bf.write_text('ExecuteAFile('+json.dumps(str((folder/'fit.bf').resolve()))+');\nLFCompute('+lfs[0]+',LF_START_COMPUTE);\nLFCompute('+lfs[0]+',recomputed);\nLFCompute('+lfs[0]+',LF_DONE_COMPUTE);\nfprintf(stdout,"RECOMPUTED_LOGL=",Format(recomputed,0,15),"\\n");\n')
        with (check/'reevaluate.log').open('w') as stream:subprocess.run([str(exe.resolve()),'CPU=1',str(bf.resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True)
        matches=re.findall(r'^RECOMPUTED_LOGL=(\S+)',(check/'reevaluate.log').read_text(),re.M)
        if len(matches)!=1:raise ValueError('Missing reevaluated likelihood')
        reevaluated=float(matches[0]);delta=reevaluated-ll
        if not math.isfinite(reevaluated):raise ValueError('Nonfinite reevaluated likelihood')
        flags=[]
        if abs(delta)>1e-6:flags.append('saved_likelihood_differs')
        if abs(saved_omega-omega)>1e-10:flags.append('saved_omega_differs')
        if not ci['LB']-1e-10<=omega<=ci['UB']+1e-10:flags.append('omega_outside_profile_interval')
        if largest_component_error>1e-8:flags.append('branch_components_not_additive')
        log=(folder/'fit.log').read_text();warnings=[line.strip() for line in log.splitlines() if 'warning' in line.lower()]
        summaries.append({'case_id':case,'translation_table':rc['translation_table'],'taxa':len(taxa),'codons':ncols,'branches':len(attr),'reported_log_likelihood':ll,'reevaluated_log_likelihood':reevaluated,'likelihood_difference':delta,'reported_omega':omega,'saved_omega':saved_omega,'omega_profile_lower':ci['LB'],'omega_profile_upper':ci['UB'],'maximum_component_additivity_error':largest_component_error,'numerical_review_flags':';'.join(flags),'warning_lines':len(warnings),'marker_copy_caveat':rc['marker_copy_caveat'],'selection_eligibility':'not_established'})
        sources.append({'case_id':case,'source_receipt_sha256':sha(rp),'reevaluation_script_sha256':sha(bf),'reevaluation_log_sha256':sha(check/'reevaluate.log')})
        if count%50==0:print('Audited',count,'of',len(paths),flush=True)
    for name,rows in [('case_audit.tsv',summaries),('branch_components.tsv',branches)]:
        with (args.output/name).open('w',newline='') as stream:
            writer=csv.DictWriter(stream,list(rows[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
    flagged=sum(bool(r['numerical_review_flags']) for r in summaries)
    result={'status':'complete_saved_fit_readback_with_review_flags' if flagged else 'passed_saved_fit_readback','scope':'completed_case_snapshot' if len(paths)<config['cases'] else 'all_executed_cases','audited_cases':len(paths),'planned_cases':config['cases'],'branches':len(branches),'cases_with_numerical_review_flags':flagged,'maximum_absolute_likelihood_difference':max(abs(r['likelihood_difference']) for r in summaries),'maximum_component_additivity_error':max(r['maximum_component_additivity_error'] for r in summaries),'source_config_sha256':sha(cp),'source_cases':sources,'script_sha256':sha(Path(__file__)),'artifacts':{name:sha(args.output/name) for name in ['case_audit.tsv','branch_components.tsv']},'interpretation':'All saved likelihoods reloaded and evaluated without optimization; dimensions, column coverage, topology, branch identities and numerical report consistency checked. This does not establish an optimum, profile-interval calibration, model adequacy, conventional dS normalization, alignment/orthology/recombination eligibility or selection.'}
    (args.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['source_cases','artifacts']},indent=2))


if __name__=='__main__':main()
