#!/usr/bin/env python3
"""Evaluate longest-dS branch likelihood slices with nuisance parameters fixed."""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import json
import math
from pathlib import Path
import re
import subprocess
from Bio import Phylo, SeqIO
from audit_genus_codon_trees import sha, table


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['fits','normalized','tree-review','plan','install','output']:parser.add_argument('--'+name,type=Path,required=True)
    a=parser.parse_args()
    if a.output.exists():raise FileExistsError('Use a new immutable output directory')
    plan=json.loads(a.plan.read_text());fit=json.loads((a.fits/'receipt.json').read_text());norm=json.loads((a.normalized/'receipt.json').read_text())
    for path,key in [(a.fits/'receipt.json','fit_receipt_sha256'),(a.normalized/'receipt.json','normalization_receipt_sha256'),(a.tree_review,'tree_review_sha256')]:
        if sha(path)!=plan[key]:raise ValueError('Source differs from plan')
    for name,digest in norm['artifacts'].items():
        if sha(a.normalized/name)!=digest:raise ValueError('Changed normalization artifact')
    expected={r['case_id']:r['receipt_sha256'] for r in fit['case_receipts']};reviews={r['case_id']:r for r in table(a.tree_review)}
    branches={}
    for r in table(a.normalized/'normalized_branches.tsv'):branches.setdefault(r['case_id'],[]).append(r)
    if set(branches)!=set(expected) or set(reviews)!=set(expected) or len(expected)!=plan['cases']:raise ValueError('Case grid differs')
    config=json.loads((a.fits/'config.json').read_text());exe=a.install/'bin/hyphy'
    if sha(exe)!=config['executable_sha256']:raise ValueError('Changed executable')
    a.output.mkdir(parents=True)
    def run(case):
        folder=a.fits/case;receipt=json.loads((folder/'receipt.json').read_text());rc=json.loads((folder/'config.json').read_text())
        if sha(folder/'receipt.json')!=expected[case] or receipt['config_sha256']!=sha(folder/'config.json'):raise ValueError('Changed fit provenance')
        for name,digest in receipt['artifacts'].items():
            if sha(folder/name)!=digest:raise ValueError('Changed saved fit artifact')
        target=sorted(branches[case],key=lambda r:(-float(r['dS_upstream_equal_alternative_convention']),r['node']))[0]
        name=target['node'];ds=float(target['dS_upstream_equal_alternative_convention'])
        tree=Phylo.read(folder/'tree.nwk','newick');nodes={n.name:n for n in tree.find_clades()};node=nodes[name];taxa={n.name for n in tree.get_terminals()}
        side={n.name for n in node.get_terminals()};split=min((tuple(sorted(side)),tuple(sorted(taxa-side))),key=lambda s:(len(s),s))
        command=rc['command'];alignment=Path(command[command.index('--alignment')+1])
        if sha(alignment)!=rc['alignment_sha256']:raise ValueError('Changed codon alignment')
        records=list(SeqIO.parse(alignment,'fasta'));coverage={};lengths=set()
        for r in records:
            seq=str(r.seq);codons=[seq[i:i+3] for i in range(0,len(seq),3)];lengths.add(len(codons))
            if any(c!='???' and (len(c)!=3 or set(c)-set('ACGT')) for c in codons):raise ValueError('Unexpected codon state')
            coverage[r.id]=sum(c!='???' for c in codons)/len(codons)
        if len(lengths)!=1 or set(coverage)!=taxa:raise ValueError('Alignment grid differs')
        saved=(folder/'fit.bf').read_text();lf=re.findall(r'^LikelihoodFunction\s+(\S+)\s*=',saved,re.M)
        pattern=r'^([A-Za-z0-9_.]+\.tree_0\.'+re.escape(name)+r'\.t)=([^;]+);'
        variables=re.findall(pattern,saved,re.M)
        if len(lf)!=1 or len(variables)!=1:raise ValueError('Ambiguous branch variable')
        variable,original_text=variables[0];original=float(original_text)
        if not math.isfinite(original) or original<=0:raise ValueError('Nonpositive target branch parameter')
        output=a.output/case;output.mkdir()
        hbl='ExecuteAFile('+json.dumps(str((folder/'fit.bf').resolve()))+');\noriginal_t='+variable+';\nLFCompute('+lf[0]+',LF_START_COMPUTE);\n'
        for factor in plan['branch_length_multipliers']:
            hbl+=variable+'=original_t*'+str(factor)+';\nLFCompute('+lf[0]+',slice_ll);\nfprintf(stdout,"SLICE\\t'+str(factor)+'\\t",Format('+variable+',0,16),"\\t",Format(slice_ll,0,16),"\\n");\n'
        hbl+=variable+'=original_t;\nLFCompute('+lf[0]+',restored_ll);\nLFCompute('+lf[0]+',LF_DONE_COMPUTE);\nfprintf(stdout,"RESTORED=",Format(restored_ll,0,16),"\\n");\n'
        (output/'slice.bf').write_text(hbl)
        with (output/'slice.log').open('w') as stream:subprocess.run([str(exe.resolve()),'CPU=1',str((output/'slice.bf').resolve())],stdout=stream,stderr=subprocess.STDOUT,check=True)
        log=(output/'slice.log').read_text();matches=re.findall(r'^SLICE\t(\S+)\t(\S+)\t(\S+)',log,re.M)
        restored=re.findall(r'^RESTORED=(\S+)',log,re.M)
        if len(matches)!=len(plan['branch_length_multipliers']) or len(restored)!=1:raise ValueError('Incomplete slice output')
        ll=float(receipt['log_likelihood'])
        if abs(float(restored[0])-ll)>1e-6:raise ValueError('Restored likelihood differs')
        rows=[]
        for expected_factor,(factor,value,evaluated) in zip(plan['branch_length_multipliers'],matches):
            factor,value,evaluated=map(float,(factor,value,evaluated))
            if factor!=expected_factor or not math.isfinite(evaluated) or not math.isclose(value,original*factor,rel_tol=1e-10,abs_tol=1e-12):raise ValueError('Wrong slice parameter or likelihood')
            if factor==1 and abs(evaluated-ll)>1e-6:raise ValueError('Baseline likelihood differs')
            rows.append({'case_id':case,'node':name,'branch_multiplier':factor,'branch_parameter':value,'normalized_dS':ds*factor,'log_likelihood':evaluated,'delta_log_likelihood_from_saved_fit':evaluated-ll})
        byfactor={r['branch_multiplier']:r for r in rows}
        summary={'case_id':case,'target_node':name,'target_is_terminal':node.is_terminal(),'split_smaller_side_taxa':';'.join(split),'target_dS':ds,'target_dN':float(target['dN_upstream_equal_alternative_convention']),'target_branch_parameter':original,'taxa':len(taxa),'codons':next(iter(lengths)),'minimum_taxon_codon_coverage':min(coverage.values()),'target_terminal_codon_coverage':coverage[name] if node.is_terminal() else '', 'delta_ll_at_half':byfactor[.5]['delta_log_likelihood_from_saved_fit'],'delta_ll_at_double':byfactor[2]['delta_log_likelihood_from_saved_fit'],'delta_ll_at_tenfold':byfactor[10]['delta_log_likelihood_from_saved_fit'],'maximum_slice_improvement':max(r['delta_log_likelihood_from_saved_fit'] for r in rows),'marker_copy_caveat':rc['marker_copy_caveat']}
        for key in ['saturated_pairwise_distance_warning_lines','parameter_boundary_warning_lines','nni_convergence_warning_lines','other_warning_lines']:summary[key]=reviews[case][key]
        proof={'case_id':case,'source_receipt_sha256':sha(folder/'receipt.json'),'slice_script_sha256':sha(output/'slice.bf'),'slice_log_sha256':sha(output/'slice.log')}
        return summary,rows,proof
    summaries=[];rows=[];proofs=[]
    with ThreadPoolExecutor(max_workers=plan['workers']) as pool:
        futures=[pool.submit(run,case) for case in sorted(expected)]
        for i,future in enumerate(as_completed(futures),1):
            summary,points,proof=future.result();summaries.append(summary);rows.extend(points);proofs.append(proof)
            if i%100==0:print('Sliced',i,flush=True)
    summaries.sort(key=lambda r:(-r['target_dS'],r['case_id']));rows.sort(key=lambda r:(r['case_id'],r['branch_multiplier']))
    for name,data in [('case_review.tsv',summaries),('likelihood_slices.tsv',rows)]:
        with (a.output/name).open('w',newline='') as stream:
            w=csv.DictWriter(stream,list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result={'status':'complete_full_case_conditional_longest_branch_slices','cases':len(summaries),'slice_evaluations':len(rows),'restored_baseline_checks':len(proofs),'plan_sha256':sha(a.plan),'script_sha256':sha(Path(__file__)),'source_fit_receipt_sha256':sha(a.fits/'receipt.json'),'normalization_receipt_sha256':sha(a.normalized/'receipt.json'),'checks':sorted(proofs,key=lambda r:r['case_id']),'artifacts':{name:sha(a.output/name) for name in ['case_review.tsv','likelihood_slices.tsv']},'interpretation':plan['interpretation']}
    (a.output/'receipt.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['checks','artifacts']},indent=2))


if __name__=='__main__':main()
