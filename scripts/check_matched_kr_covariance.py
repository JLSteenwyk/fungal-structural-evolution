"""Compare adjustment against pinned author R code at identical input values."""
import itertools
import json
from pathlib import Path
import subprocess
import numpy as np
from threadpoolctl import threadpool_limits
from ancestral_chain_attempt import sha, write_json
from matched_covariance_information_fast import covariance_information
from matched_kr_covariance import adjust_covariance


def main():
    threadpool_limits(1)
    provenance = Path('metadata/pbkrtest_covariance_reference_source_20260928.json')
    source = json.loads(provenance.read_text())
    assert sha(source['path']) == source['sha256']
    for path, rec in source['additional_files'].items():
        assert sha(path) == rec['sha256']
    out = Path('results/model_validation/matched-kr-covariance-checks-20260928-v1')
    out.mkdir(parents=True, exist_ok=False)
    rng = np.random.default_rng(20260928)
    n = 48; bg = np.repeat(np.arange(12), 4); fam = bg // 3
    x = np.column_stack([np.ones(n), rng.normal(size=(n, 3))])
    f = rng.normal(size=(n, 5))/3
    kernels = [np.eye(n), (bg[:, None] == bg).astype(float),
               (fam[:, None] == fam).astype(float), f@f.T]
    np.savetxt(out/'X.csv', x, delimiter=',')
    for i, k in enumerate(kernels):
        np.savetxt(out/f'G{i}.csv', k, delimiter=',')
    results = []
    for ratios in itertools.product([0., .7], repeat=3):
        for scale in [.2, 3.]:
            v = scale*np.array([1., *ratios])
            c = covariance_information(bg, fam, f, x, v)
            adjusted = adjust_covariance(c)
            assert adjusted['status'] == 'candidate_adjustment_pending_statistical_validation'
            idx = len(results)
            np.savetxt(out/f'V{idx}.csv', sum(a*k for a,k in zip(v,kernels)), delimiter=',')
            np.savetxt(out/f'Phi{idx}.csv', c['conditional_beta_covariance'], delimiter=',')
            results.append(adjusted)
    # Run unmodified author covariance routine and its original index helper.
    script = '''library(Matrix)
library(MASS)
args <- commandArgs(trailingOnly=TRUE)
source(args[1]); source(args[2]); root <- args[3]
read <- function(name) as.matrix(read.csv(file.path(root,name), header=FALSE))
X <- read("X.csv")
G <- lapply(0:3, function(i) Matrix(read(paste0("G",i,".csv"))))
for (i in 0:15) {
  z <- vcovAdj_internal(Matrix(read(paste0("Phi",i,".csv"))),
      list(Sigma=Matrix(read(paste0("V",i,".csv"))), G=G, n.ggamma=4), X)
  write.table(as.matrix(z), file.path(root,paste0("adjusted",i,".csv")), sep=",", row.names=FALSE, col.names=FALSE)
  write.table(as.matrix(attr(z,"W")), file.path(root,paste0("W",i,".csv")), sep=",", row.names=FALSE, col.names=FALSE)
}
writeLines(c(R.version.string, paste("Matrix",packageVersion("Matrix")), paste("MASS",packageVersion("MASS"))), file.path(root,"R_versions.txt"))
'''
    (out/'reference.R').write_text(script)
    helper = str(Path(source['path']).parent/'KR_utils.R')
    subprocess.run(['Rscript', str(out/'reference.R'), helper, source['path'], str(out)], check=True)
    errors = []
    for i, result in enumerate(results):
        for key, prefix in [('adjusted_beta_covariance','adjusted'), ('variance_parameter_covariance','W')]:
            expected = np.loadtxt(out/f'{prefix}{i}.csv', delimiter=',')
            np.testing.assert_allclose(result[key], expected, rtol=1e-9, atol=1e-10)
            errors.append(float(np.max(abs(result[key]-expected))))
    confounded = (fam[:, None] == np.unique(fam)).astype(float)
    bad = covariance_information(bg, fam, confounded, x, [1., .2, .3, .4])
    assert adjust_covariance(bad)['status'] == 'information_requires_review'
    record = dict(status='passed_covariance_adjustment_author_code_comparison', cases=16,
        maximum_absolute_disagreement=max(errors), singular_information_retained_for_review=True,
        source_receipt_sha256=sha(provenance), source_commit=source['commit'],
        pins={str(p):sha(p) for p in [Path(__file__),Path('scripts/matched_kr_covariance.py'),
             Path('scripts/matched_covariance_information_fast.py'),Path('scripts/matched_covariance_information.py'),
             Path('scripts/matched_mixed_covariance.py')]},
        artifacts={p.name:sha(p) for p in out.iterdir()},
        scope='Identical fixed variance/design inputs; unmodified author covariance routine, '
              'not installed pbkrtest or an independent fit. Formula agreement including zero '
              'components does not validate boundary coverage, degrees of freedom or intervals.')
    write_json(out/'receipt.json', record)
    write_json(Path('metadata/matched_kr_covariance_checks_20260928.json'), record)
    print(json.dumps({k:v for k,v in record.items() if k not in ['pins','artifacts']}, indent=2))


if __name__ == '__main__':
    main()
