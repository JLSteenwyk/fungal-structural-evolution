#!/usr/bin/env python3
"""Change only the starting alignment of immutable BAli-Phy model programs.

The native distribution, its annotated density, modifiable representation and
transition kernels remain the installed originals. Never fixes ancestral
lengths, observes the reference matrix, clips parameters or changes a prior.
"""
import collections
import hashlib
import json
import re
from pathlib import Path

from Bio import SeqIO


HELPER = r'''
-- BEGIN project reference alignment initialization v1
referenceAlignmentStart tree branchHMMs tipLengths referenceData = do
  AlignmentOnTree _ nodeCount fixedTipLengths _ <- sample_alignment tree branchHMMs tipLengths
  let AlignmentOnTree _ _ _ referenceBranches = alignmentOnTreeFromSequences tree referenceData
  return (branchHMMs `deepseq` (AlignmentOnTree tree nodeCount fixedTipLengths referenceBranches))

sampleReferenceAlignmentWithProps dist@(PhyloAlignment tree imodel tipLengths branchHMMs) referenceData = do
  x <- RanDistribution3 dist alignment_effect triggeredModifiableAlignment
         (referenceAlignmentStart tree branchHMMs tipLengths referenceData)
  let props = getProperties' x dist
      getProperties' :: a -> d -> DistProperties d
      getProperties' x _ = Effect.getProperties x
  return (x, props)

referenceInitializationAudit alignment sequenceData properties =
  ["referenceInitialization" %>%
    ["extantAlignment" %=% (toFasta (align alignment sequenceData)),
     "nodeCount" %=% nodeCount,
     "fixedLengthNodes" %=% (IntMap.keys fixedTipLengths),
     "allLengthNodes" %=% (IntMap.keys (mkSequenceLengthsMap alignment)),
     "priorFromDistribution" %=% (ln (probability properties)),
     "priorFromDirectFormula" %=% directPrior]]
  where
    AlignmentOnTree tree nodeCount fixedTipLengths branches = alignment
    directPrior = ln (alignment_pr_top branches tree (hmms properties)) +
      sum [(1 - fromIntegral (nodeDegree tree node)) * ln (lengthp properties (sequenceLength alignment node))
           | node <- internalNodes tree]
-- END project reference alignment initialization v1
'''
IMPORTS = ('import Probability.Distribution.PhyloAlignment\n'
           'import Control.DeepSeq\nimport qualified Effect\n')


def substitutions(original):
    loads = re.findall(r';sequenceData <- \(mkUnalignedCharacterData aa\) <\$> \(loadSequences ("[^"\n]+")\)', original)
    assert len(loads) == 1, 'Expected exactly one immutable aligned input load'
    return [
        ('import Probability.Random\n', 'import Probability.Random\n' + IMPORTS),
        ('\nmodel sequenceData tree ', HELPER + '\nmodel sequenceData referenceData initializationAudit tree '),
        (';(alignment,properties_A) <- sampleWithProps (phyloAlignment tree imodel scale1 sequence_lengths)',
         ';(alignment,properties_A) <- sampleReferenceAlignmentWithProps (phyloAlignment tree imodel scale1 sequence_lengths) referenceData'),
        (';return (parameterLogValues loggerValues)',
         ';return (parameterLogValues loggerValues ++ (if initializationAudit then referenceInitializationAudit alignment sequenceData properties_A else []))'),
        (';sequenceData <- (mkUnalignedCharacterData aa) <$> (loadSequences ' + loads[0] + ')',
         ';inputSequences <- loadSequences ' + loads[0] + '\n'
         ';let {sequenceData = mkUnalignedCharacterData aa inputSequences}\n'
         ';let {referenceData = mkAlignedCharacterData aa inputSequences}'),
        ('(model sequenceData tree logParamsTSV', '(model sequenceData referenceData isTest tree logParamsTSV'),
    ]


def transform(original):
    """Require each expected edit exactly once and prove byte-exact reversal."""
    edits = substitutions(original)
    result = original
    for before, after in edits:
        assert result.count(before) == 1, before
        assert after not in result, after
        result = result.replace(before, after)
    restored = result
    for before, after in reversed(edits):
        assert restored.count(after) == 1, after
        restored = restored.replace(after, before)
    assert restored == original
    assert result.count('observe sequenceData (phyloCTMC tree alignment smodel scale1)') == 1
    assert 'observe referenceData' not in result
    return result


def alignment_records(path):
    records = [(r.id, str(r.seq).upper()) for r in SeqIO.parse(str(path), 'fasta')]
    assert records and len({r[0] for r in records}) == len(records), 'Missing or duplicate labels'
    assert len({len(r[1]) for r in records}) == 1, 'Unequal aligned sequence lengths'
    assert all(set(s) <= set('ARNDCQEGHILKMFPSTWYVBZX?-') for _, s in records), 'Unexpected residue symbol'
    return records


def homology_signature(records):
    """Canonical column-to-tip-residue memberships; ignore all-gap columns/order.

    Native tree projection may reorder independent insertion columns. This
    compares homology and residue identities without demanding identical text.
    """
    assert len({label for label, _ in records}) == len(records)
    assert len({len(s) for _, s in records}) == 1
    positions = collections.Counter()
    columns = []
    for offset in range(len(records[0][1])):
        members = []
        for label, sequence in records:
            residue = sequence[offset]
            if residue not in '-':
                positions[label] += 1
                members.append((label, positions[label], residue))
        if members:
            columns.append(tuple(sorted(members)))
    return sorted(columns)


def validate_initial_json(value, reference, tips):
    """Native initial audit is conditional model/startup proof, not convergence."""
    import io
    import math
    audit = value['parameters/']['referenceInitialization/']
    records = [(x.id, str(x.seq).upper()) for x in SeqIO.parse(io.StringIO(audit['extantAlignment']), 'fasta')]
    assert {x[0] for x in records} == set(tips) == {x[0] for x in reference}
    assert homology_signature(records) == homology_signature(reference), 'Reference homology changed'
    fixed, all_nodes = audit['fixedLengthNodes'], audit['allLengthNodes']
    assert len(fixed) == len(set(fixed)) == len(tips)
    assert len(all_nodes) == len(set(all_nodes)) == audit['nodeCount']
    assert set(fixed) < set(all_nodes), 'Ancestral lengths were fixed'
    scores = [value[name] for name in ['prior', 'likelihood', 'posterior']]
    assert all(math.isfinite(x) for x in scores)
    assert math.isclose(scores[0] + scores[1], scores[2], abs_tol=1e-7, rel_tol=1e-12)
    a, b = audit['priorFromDistribution'], audit['priorFromDirectFormula']
    assert math.isfinite(a) and math.isfinite(b)
    assert math.isclose(a, b, abs_tol=1e-8, rel_tol=1e-12), 'Native alignment densities disagree'
    return dict(tips=len(tips), nodes=audit['nodeCount'], fixed_tip_lengths=len(fixed),
                ancestral_lengths_free=len(all_nodes) - len(fixed),
                extant_homology_columns=len(homology_signature(reference)),
                prior=scores[0], likelihood=scores[1], posterior=scores[2],
                alignment_log_prior=a, direct_alignment_log_prior=b,
                reference_homology_preserved=True, scientific_eligibility=False)


def digest(value):
    return hashlib.sha256(value.encode()).hexdigest()
