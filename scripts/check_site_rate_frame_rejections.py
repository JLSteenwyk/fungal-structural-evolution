#!/usr/bin/env python3
"""Check frame readback rejects rehashed semantic corruptions of a completed cohort."""
import csv,json,subprocess,sys,tempfile
from pathlib import Path
from audit_busco_gene_copies import sha


def main():
    frame=Path('results/phylogeny/site-rate-exposure-frame-esmfold-all-completed-20260923-v1')
    base=[sys.executable,'scripts/readback_site_rate_exposure_frame.py','--rates','results/phylogeny/paired-rate-heterogeneity-optimized-esmfold-all-completed-20260923-v1','--rate-readback','results/phylogeny/paired-rate-heterogeneity-readback-esmfold-all-completed-20260923-v1/receipt.json','--exposure','results/phylogeny/site-parsimony-exposure-esmfold-all-completed-20260923-v1','--exposure-readback','results/phylogeny/site-parsimony-exposure-esmfold-all-completed-20260923-v1-readback/receipt.json','--inputs','results/phylogeny/paired-inputs-esmfold-all-completed-20260922-v1']
    with (frame/'site_rate_exposure.tsv').open() as f:reader=csv.DictReader(f,delimiter='\t');fields=reader.fieldnames;original=list(reader)
    rejected=[]
    with tempfile.TemporaryDirectory() as temp:
        for name,message in [('rate','Rate value differs'),('entropy','Derived covariate differs'),('missing','Frame universe differs')]:
            root=Path(temp)/name;root.mkdir();rows=[dict(x) for x in original]
            if name=='rate':rows[0]['aa_gamma_rate']=str(float(rows[0]['aa_gamma_rate'])+1)
            elif name=='entropy':rows[0]['aa_entropy_nats']=str(float(rows[0]['aa_entropy_nats'])+.1)
            else:rows.pop()
            with (root/'site_rate_exposure.tsv').open('w') as f:
                writer=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(rows)
            (root/'fit_diagnostics.tsv').symlink_to((frame/'fit_diagnostics.tsv').resolve())
            receipt=json.loads((frame/'receipt.json').read_text());receipt['artifacts']={n:sha(root/n) for n in receipt['artifacts']};(root/'receipt.json').write_text(json.dumps(receipt))
            result=subprocess.run(base+['--frame',str(root),'--output',str(root/'checked.json')],capture_output=True,text=True)
            assert result.returncode!=0 and message in result.stderr and not (root/'checked.json').exists(),(name,result.stderr)
            rejected.append(name)
    record=dict(status='passed_site_rate_frame_semantic_rejections',rejected_cases=rejected,source_frame_receipt_sha256=sha(frame/'receipt.json'),checker_sha256=sha(Path('scripts/readback_site_rate_exposure_frame.py')),test_sha256=sha(Path(__file__)),scope='Temporary copies of the completed cohort with recomputed artifact hashes; rejects altered rates, altered entropy and missing site. Original production data unchanged.')
    Path('metadata/site_rate_frame_rejection_checks_20260927.json').write_text(json.dumps(record,indent=2)+'\n');print(json.dumps(record,indent=2))


if __name__=='__main__':main()
