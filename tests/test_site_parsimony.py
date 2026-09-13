import io
import itertools
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from Bio import Phylo
from summarize_site_parsimony_exposure import minimum_changes


def exhaustive(tree,observed):
    nodes=list(tree.find_clades());free=[n for n in nodes if not n.is_terminal() or observed[n.name]=='?']
    fixed={n:observed[n.name] for n in nodes if n not in free};best=100
    for assignment in itertools.product('AB',repeat=len(free)):
        states=fixed|dict(zip(free,assignment))
        best=min(best,sum(states[n]!=states[c] for n in nodes for c in n.clades))
    return best


class ParsimonyTest(unittest.TestCase):
    def test_exhaustive_missing_and_multifurcating_trees(self):
        for newick in ['((a,b),(c,d));','(a,b,c,d);','(a,(b,(c,d)));']:
            tree=Phylo.read(io.StringIO(newick),'newick')
            for pattern in itertools.product('AB?',repeat=4):
                observed=dict(zip('abcd',pattern))
                with self.subTest(tree=newick,pattern=pattern):
                    self.assertEqual(minimum_changes(tree,observed,'AB').tolist(),[exhaustive(tree,observed)])
    def test_hard_polytomy_not_silently_resolved(self):
        tree=Phylo.read(io.StringIO('(a,b,c,d);'),'newick')
        self.assertEqual(minimum_changes(tree,dict(zip('abcd','AABB')),'AB').tolist(),[2])
    def test_distinct_and_unknown_columns(self):
        tree=Phylo.read(io.StringIO('(a,b,c);'),'newick')
        self.assertEqual(minimum_changes(tree,{'a':'A?A','b':'B?A','c':'C?A'},'ABC').tolist(),[2,0,0])
    def test_invalid_taxa_or_symbols(self):
        tree=Phylo.read(io.StringIO('(a,b);'),'newick')
        for seqs in [{'a':'A'},{'a':'A','b':'X'}]:
            with self.assertRaises(ValueError):minimum_changes(tree,seqs,'AB')

if __name__=='__main__':unittest.main()
