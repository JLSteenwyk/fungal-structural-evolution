#!/usr/bin/env python3
"""Prepare isolated precision-only and MP-scaling cache-refresh sensitivity builds."""
import argparse
import difflib
import json
from pathlib import Path
import shutil
from ancestral_chain_attempt import sha, write_json


def replace_once(text, old, new):
    assert text.count(old) == 1, old
    return text.replace(old, new)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output', type=Path, required=True)
    args = ap.parse_args()
    source = Path('data/software_audits/fastml-3.11/source/FastML.v3.11')
    assert not args.output.exists()
    inventory = {str(p.relative_to(source)): sha(p) for p in source.rglob('*') if p.is_file()}
    variants = []
    args.output.mkdir(parents=True)
    for label in ['precision-only', 'precision-cache-refresh']:
        target = args.output / label
        shutil.copytree(source, target)
        changes = {}
        for name in ['gainLoss.cpp', 'gainLossUtils.cpp']:
            relative = Path('programs/gainLoss') / name
            old = (source / relative).read_text()
            text = old
            if name == 'gainLoss.cpp':
                text = replace_once(text, 'ofstream modelParamsStream(modelParams.c_str());',
                                    'ofstream modelParamsStream(modelParams.c_str());\n\tmodelParamsStream.precision(17);')
                text = replace_once(text, 'AncestralReonstructPosteriorStream.precision(PRECISION);',
                                    'AncestralReonstructPosteriorStream.precision(17);')
                if label == 'precision-cache-refresh':
                    start = text.index('void gainLoss::multipleAllBranchesByFactorAtStartByMaxParsimonyCost(')
                    end = text.index('\n}', start) + 2
                    function = text[start:end]
                    function = replace_once(function, '\t_tr.multipleAllBranchesByFactor(factorBL);',
                        '\t_tr.multipleAllBranchesByFactor(factorBL);\n'
                        '\t// Project sensitivity: refresh ascertainment after MP branch scaling.\n'
                        '\tif(_unObservableData_p){\n'
                        '\t\tif(gainLossOptions::_gainLossDist)\n'
                        '\t\t\t_unObservableData_p->setLforMissingData(_tr,_spVVec,_gainDist,_lossDist);\n'
                        '\t\telse\n'
                        '\t\t\t_unObservableData_p->setLforMissingData(_tr,_sp);\n'
                        '\t}\n'
                        '\tprintTreeLikelihoodAllPosAlphTheSame(); // refresh optimizer baseline\n')
                    text = text[:start] + function + text[end:]
            else:
                text = replace_once(text, 'void printTree (tree &tr,ostream &out){',
                                    'void printTree (tree &tr,ostream &out){\n\tout.precision(17);')
                text = replace_once(text, 'ofstream treeStream(gainLossOptions::_treeOutFile.c_str());',
                                    'ofstream treeStream(gainLossOptions::_treeOutFile.c_str());\n\ttreeStream.precision(17);')
            # Preserve CRLF source convention, keeping the scientific diff small.
            (target / relative).write_bytes(text.replace('\n', '\r\n').encode())
            patch = args.output / (label + '-' + name + '.patch')
            patch.write_text(''.join(difflib.unified_diff(old.splitlines(True), text.splitlines(True),
                             fromfile='original/' + str(relative), tofile=label + '/' + str(relative))))
            changes[str(relative)] = dict(original_sha256=inventory[str(relative)],
                modified_sha256=sha(target / relative), patch=str(patch), patch_sha256=sha(patch))
        variants.append(dict(label=label, directory=str(target), changes=changes,
            build_commands=[['make', '-C', str(target / 'programs/gainLoss'), 'clean'],
                            ['make', '-C', str(target / 'programs/gainLoss'), '-j2']]))
    assert all(sha(source / path) == digest for path, digest in inventory.items())
    write_json(args.output / 'preparation.json', dict(status='prepared_not_compiled_or_qualified',
        source=str(source), source_inventory=inventory, variants=variants, script_sha256=sha(__file__),
        resources=dict(cpus=2, memory_gib=8, output_allowance_gib=1, planning_hours=[0.1, 2],
            basis='Two gainLoss rebuilds from frozen source; retain and bind the existing phylogeny static library.',
            gpu=False, paid_resources=False),
        scope='Precision-only changes serialization. Cache-refresh additionally refreshes exclusion '
              'probability and optimizer baseline after MP initial scaling. Other optimizer/cache '
              'routes are not certified; scientific qualification requires paired runs and replay.'))


if __name__ == '__main__':
    main()
