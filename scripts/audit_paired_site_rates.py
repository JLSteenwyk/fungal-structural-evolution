#!/usr/bin/env python3
"""Audit fixed-topology rate exports, category bounds and site-likelihood accounting."""
import argparse,json,math,re
from decimal import Decimal
from pathlib import Path
from Bio import SeqIO
from assess_pae_sensitivity import checked_receipt
from audit_busco_gene_copies import sha,read_table
from prepare_paired_phylogenetic_inputs import write_table
from run_paired_marker_fits import tree_edges


def rounding_bound(text):return .5*10.**Decimal(text).as_tuple().exponent


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['rates','inputs','fits','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if a.output.exists():raise FileExistsError('Use immutable output')
    r=json.loads((a.rates/'receipt.json').read_text());c=json.loads((a.rates/'config.json').read_text());checked_receipt(a.inputs)
    if r['config_sha256']!=sha(a.rates/'config.json') or any(c['source_receipts'][k]!=sha(getattr(a,k)/'receipt.json') for k in ['inputs','fits']):raise ValueError('Source lineage differs')
    for path,h in c['pinned_files'].items():
        if sha(Path(path))!=h:raise ValueError('Pinned source changed')
    ready=[x for x in read_table(a.inputs/'marker_summary.tsv') if x['status']=='ready_for_inference'];labels=['aa','3di_af','3di_af_empirical','3di_llm']
    expected={x['marker']+'/'+label for x in ready for label in labels}
    if set(r['fit_receipts'])!=expected:raise ValueError('Incomplete fit grid')
    fits=[];sites=[];warnings=[]
    for m in ready:
        marker=m['marker'];width=int(m['retained_columns']);taxa={x.id for x in SeqIO.parse(a.inputs/marker/'aa.faa','fasta')};topology=a.fits/marker/'aa.treefile'
        for label in labels:
            folder=a.rates/marker;rp=folder/(label+'.receipt.json');fr=json.loads(rp.read_text());request=json.loads((folder/(label+'.config.json')).read_text())
            if sha(rp)!=r['fit_receipts'][marker+'/'+label] or fr['parent_config_sha256']!=r['config_sha256'] or request['parent_config_sha256']!=r['config_sha256'] or request['topology_sha256']!=sha(topology):raise ValueError('Fit provenance differs')
            for name,h in fr['artifacts'].items():
                if sha(folder/name)!=h:raise ValueError('Fit artifact changed')
            alignment=a.inputs/marker/('aa.faa' if label=='aa' else '3di.faa')
            if request['alignment_sha256']!=sha(alignment):raise ValueError('Alignment changed')
            report=(folder/(label+'.iqtree')).read_text();model=request['command'][request['command'].index('-m')+1]
            heterogeneity=c.get('heterogeneity','G4')
            original_command=json.loads((a.fits/marker/(label+'.config.json')).read_text())['command']
            original_model=original_command[original_command.index('-m')+1]
            if heterogeneity not in ['G4','R4'] or not original_model.endswith('+G4') or model!=original_model[:-2]+heterogeneity:raise ValueError('Unexpected model transformation')
            rate_name={'G4':'Gamma','R4':'FreeRate'}[heterogeneity]
            if 'Model of substitution: '+model+'\n' not in report or 'Model of rate heterogeneity: '+rate_name+' with 4 categories' not in report:raise ValueError('Model identity differs')
            if set(tree_edges(folder/(label+'.treefile'),taxa))!=set(tree_edges(topology,taxa)):raise ValueError('Topology changed')
            block=report.split('Category  Relative_rate  Proportion',1)[1].strip().split('\n\n',1)[0].split('Relative rates',1)[0]
            category={int(x[0]):(float(x[1]),float(x[2])) for line in block.splitlines() if len(x:=line.split())==3 and x[0].isdigit()}
            if set(category)!={1,2,3,4} or any(not math.isfinite(v) or v<0 for pair in category.values() for v in pair) or not math.isclose(sum(x[1] for x in category.values()),1,abs_tol=.00021):raise ValueError('Invalid Gamma categories')
            raw=[x.split() for x in (folder/(label+'.rate')).read_text().splitlines() if x.strip() and not x.startswith('#')]
            header=raw[0];data=raw[1:]
            if header!=['Site','Rate','Cat','C_Rate'] or len(data)!=width or [int(x[0]) for x in data]!=list(range(1,width+1)):raise ValueError('Rate output grid differs')
            for entry in data:
                index,rate,cat,catrate=int(entry[0]),float(entry[1]),int(entry[2]),float(entry[3])
                if any(not math.isfinite(x) or x<0 for x in [rate,catrate]) or cat not in category or abs(catrate-category[cat][0])>.0001 or not min(x[0] for x in category.values())-.0001<=rate<=max(x[0] for x in category.values())+.0001:raise ValueError('Site rate/category bound failed')
                sites.append({'marker':marker,'fit':label,'paired_column_1based':index,'posterior_mean_relative_rate':rate,'modal_rate_category':cat,'modal_category_relative_rate':catrate})
            lh=(folder/(label+'.sitelh')).read_text().split();lltext=re.search(r'Log-likelihood of the tree: ([\d.eE+-]+)',report).group(1)
            if lh[:2]!=['1',str(width)] or len(lh)!=width+3:raise ValueError('Site likelihood grid differs')
            values=[float(x) for x in lh[3:]]
            if any(not math.isfinite(x) or x>0 for x in values):raise ValueError('Invalid site likelihood')
            error=abs(math.fsum(values)-float(lltext));bound=sum(rounding_bound(x) for x in lh[3:])+rounding_bound(lltext)+1e-8
            if error>bound:raise ValueError('Site likelihood sum differs beyond printed precision')
            oldreport=(a.fits/marker/(label+'.iqtree')).read_text();old=float(re.search(r'Log-likelihood of the tree: ([\d.eE+-]+)',oldreport).group(1))
            warning=sorted(set(line.strip() for line in (report+'\n'+(folder/(label+'.log')).read_text()).splitlines() if re.match(r'\s*(WARNING|ERROR):',line)))
            warnings.extend({'marker':marker,'fit':label,'warning':line} for line in warning)
            fits.append({'marker':marker,'fit':label,'sites':width,'log_likelihood':float(lltext),'original_fit_log_likelihood':old,'log_likelihood_change':float(lltext)-old,'site_likelihood_sum_error':error,'printed_precision_bound':bound,'gamma_alpha':float(re.search(r'Gamma shape alpha: ([\d.eE+-]+)',report).group(1)) if heterogeneity=='G4' else '', 'rate_heterogeneity':heterogeneity,'warnings':len(warning)})
    if len(fits)!=r['fits'] or len(sites)!=r['site_rate_rows']:raise ValueError('Aggregate dimensions differ')
    a.output.mkdir(parents=True)
    for name,rows in [('fit_summary.tsv',fits),('site_rates.tsv',sites),('warnings.tsv',warnings)]:
        if rows:write_table(a.output/name,rows)
    out={'status':'passed_full_site_rate_output_audit','rate_receipt_sha256':sha(a.rates/'receipt.json'),'script_sha256':sha(Path(__file__)),'markers':len(ready),'fits':len(fits),'site_rate_rows':len(sites),'fits_with_warnings':sum(x['warnings']>0 for x in fits),'scope':'All fit provenance, model/topology identities, rate grids, Gamma-category bounds and site-likelihood sums checked at exported precision. Posterior site rates themselves were not independently recomputed; these are conditional empirical-Bayes estimates, not calibrated confidence or selection tests.','artifacts':{f.name:sha(f) for f in a.output.iterdir()}}
    out['scope']=out['scope'].replace('Gamma-category','four-category rate/weight')
    (a.output/'receipt.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))

if __name__=='__main__':main()
