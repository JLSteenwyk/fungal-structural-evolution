#!/usr/bin/env python3
"""Bind every completed MG94 case to its independently audited nucleotide tree."""
import argparse
import json
from pathlib import Path
from Bio import Phylo
from audit_genus_codon_trees import sha, split_map, table


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ['fits','trees','audit','information','output']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use a new binding receipt path')
    fit=json.loads((a.fits/'receipt.json').read_text());audit=json.loads((a.audit/'receipt.json').read_text())
    if fit['status']!='complete_full_genus_mg94_execution_pending_audit' or audit['status']!='passed_full_genus_tree_audit':raise ValueError('Complete source stages required')
    config=json.loads((a.fits/'config.json').read_text())
    if fit['config_sha256']!=sha(a.fits/'config.json') or config['tree_config_sha256']!=sha(a.trees/'config.json') or audit['source_config_sha256']!=sha(a.trees/'config.json') or config['information_sha256']!=sha(a.information):raise ValueError('Changed root provenance')
    expected={r['case_id'] for r in table(a.information) if r['status']=='ready_for_supported_tree_diagnostic'}
    fit_cases={r['case_id']:r['receipt_sha256'] for r in fit['case_receipts']}
    tree_cases={r['case_id']:r['receipt_sha256'] for r in audit['source_case_receipts']}
    if set(fit_cases)!=expected or set(tree_cases)!=expected or len(fit['case_receipts'])!=len(expected) or len(audit['source_case_receipts'])!=len(expected):raise ValueError('Incomplete, duplicate or extra case binding')
    for case in sorted(expected):
        f=a.fits/case;t=a.trees/case
        if sha(f/'receipt.json')!=fit_cases[case] or sha(t/'receipt.json')!=tree_cases[case]:raise ValueError('Changed source case receipt')
        rc=json.loads((f/'config.json').read_text());fr=json.loads((f/'receipt.json').read_text());tr=json.loads((t/'receipt.json').read_text())
        if fr['config_sha256']!=sha(f/'config.json') or rc['source_tree_receipt_sha256']!=tree_cases[case] or rc['parent_config_sha256']!=sha(a.fits/'config.json'):raise ValueError('Case provenance does not match audited tree')
        if sha(t/'tree.treefile')!=tr['artifacts']['tree.treefile'] or sha(f/'tree.nwk')!=rc['tree_sha256']:raise ValueError('Changed tree bytes')
        original=Phylo.read(t/'tree.treefile','newick');exported=Phylo.read(f/'tree.nwk','newick')
        taxa={n.name for n in original.get_terminals()}
        if set(split_map(original,taxa))!=set(split_map(exported,taxa)):raise ValueError('Sanitized topology differs from audited source')
    result={'status':'passed_full_mg94_to_audited_tree_binding','cases':len(expected),'source_fit_receipt_sha256':sha(a.fits/'receipt.json'),'source_tree_audit_sha256':sha(a.audit/'receipt.json'),'information_sha256':sha(a.information),'script_sha256':sha(Path(__file__)),'scope':'Every expected information-screened case occurs exactly once in the completed MG94 grid and complete tree audit, with matching source case hashes and identical input topologies after label sanitization. Does not audit fitted likelihoods or confer biological eligibility.'}
    a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))


if __name__=='__main__':main()
