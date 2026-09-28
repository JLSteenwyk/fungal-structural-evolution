"""Compare scalar KR moments and candidate limits with unmodified author R code."""
import json
import itertools
from pathlib import Path
import subprocess
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from matched_kr_scalar import scalar_interval


def main():
    threadpool_limits(1)
    root = Path('results/model_validation/matched-kr-covariance-checks-20260928-v1')
    receipt = json.loads((root/'receipt.json').read_text())
    for name, digest in receipt['artifacts'].items():
        assert sha(root/name) == digest
    source = json.loads(Path('metadata/pbkrtest_scalar_reference_source_20260928.json').read_text())
    assert sha(source['path']) == source['sha256']
    out = Path('results/model_validation/matched-kr-scalar-checks-20260928-v2')
    out.mkdir(parents=True, exist_ok=False)
    wrapper = '''library(Matrix)
library(MASS)
a <- commandArgs(trailingOnly=TRUE)
source(file.path(a[1],"KR_utils.R")); source(file.path(a[1],"KR_vcovAdj.R"))
expr <- parse(file.path(a[1],"KR_modcomp.R"))
for (e in expr) if (is.call(e) && identical(e[[1]],as.name("<-")) && identical(e[[2]],as.name(".KR_adjust"))) eval(e)
read <- function(n) as.matrix(read.csv(file.path(a[2],n),header=FALSE))
X <- read("X.csv"); G <- lapply(0:3,function(i) Matrix(read(paste0("G",i,".csv"))))
beta <- c(.1,.2,.3,.4); contrasts <- rbind(c(1,0,0,0),c(1,.2,-.3,.4))
rows <- NULL
for(i in 0:15) {
  Phi <- Matrix(read(paste0("Phi",i,".csv")))
  z <- vcovAdj_internal(Phi,list(Sigma=Matrix(read(paste0("V",i,".csv"))),G=G,n.ggamma=4),X)
  for(j in 1:2) {
    L <- matrix(contrasts[j,],nrow=1)
    ans <- .KR_adjust(z,Phi,L,beta,rep(0,4))
    variance <- as.numeric(L %*% z %*% t(L))
    half <- qt(.975,ans$ddf)*sqrt(variance/ans$F.scaling)
    estimate <- as.numeric(L %*% beta)
    rows <- rbind(rows,c(i,j,ans$ddf,ans$F.scaling,variance,estimate-half,estimate+half))
  }
}
write.table(rows,file.path(a[3],"reference.csv"),sep=",",row.names=FALSE,col.names=FALSE)
'''
    (out/'reference.R').write_text(wrapper)
    subprocess.run(['Rscript',str(out/'reference.R'),str(Path(source['path']).parent),str(root),str(out)],check=True)
    expected = np.loadtxt(out/'reference.csv',delimiter=',')
    x = np.loadtxt(root/'X.csv',delimiter=',')
    kernels = [np.loadtxt(root/f'G{i}.csv',delimiter=',') for i in range(4)]
    contrasts = [np.array([1.,0,0,0]),np.array([1.,.2,-.3,.4])]
    results = []; errors = []
    for i in range(16):
        inverse = np.linalg.inv(np.loadtxt(root/f'V{i}.csv',delimiter=','))
        t = inverse@x; phi = np.linalg.inv(x.T@t); p = inverse-t@phi@t.T
        c = dict(conditional_beta_covariance=phi,
                 first_contractions=np.array([-t.T@g@t for g in kernels]),
                 second_contractions=np.array([[t.T@g@inverse@h@t for h in kernels] for g in kernels]),
                 expected_reml_information=np.array([[.5*np.trace(p@g@p@h) for h in kernels] for g in kernels]),
                 exact_zero_components=np.array([False, *[v == 0 for v in list(itertools.product([0., .7], repeat=3))[i//2]]]))
        for j,l in enumerate(contrasts):
            result = scalar_interval(c,[.1,.2,.3,.4],l)
            assert result['status'] == 'candidate_interval_pending_coverage_validation'
            assert result['exact_zero_components'] == c['exact_zero_components'].tolist()
            actual = [result[k] for k in ['denominator_df','f_scaling','adjusted_variance','lower','upper']]
            np.testing.assert_allclose(actual,expected[2*i+j,2:],rtol=1e-9,atol=1e-10)
            errors.append(float(np.max(abs(np.array(actual)-expected[2*i+j,2:]))))
            results.append(result)
    rejected = 0
    for contrast,alpha in [(np.zeros(4),.05),(np.ones(3),.05),(np.ones(4),0),(np.ones(4),1)]:
        try: scalar_interval(c,[.1,.2,.3,.4],contrast,alpha)
        except ValueError: rejected += 1
        else: raise AssertionError('Invalid contrast or level accepted')
    write_json(out/'candidates.json',results)
    record = dict(status='passed_32_scalar_author_code_comparisons',cases=len(results),
        maximum_absolute_disagreement=max(errors),invalid_inputs_rejected=rejected,
        source_receipt_sha256=sha('metadata/pbkrtest_scalar_reference_source_20260928.json'),
        covariance_fixture_receipt_sha256=sha(root/'receipt.json'),
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/matched_kr_scalar.py'),Path('scripts/matched_kr_covariance.py')]},
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Scalar contrasts only; same synthetic variances/designs as covariance checks, '
              'including boundary models. Dense Python contractions compared with author R '
              'moments and t limits. No refitting, coverage, joint-test or multiplicity claim.')
    write_json(out/'receipt.json',record)
    write_json(Path('metadata/matched_kr_scalar_checks_20260928.json'),record)
    print({k:v for k,v in record.items() if k not in ['pins','artifacts']})


if __name__ == '__main__':
    main()
