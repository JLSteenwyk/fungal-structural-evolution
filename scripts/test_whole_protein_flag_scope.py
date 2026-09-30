#!/usr/bin/env python3
"""Exercise complete scope, explicit errors, and rejection of inconsistent gates."""
import unittest
from whole_protein_flag_followup import full_scope, PASS, FLAG, ERROR


def rows(ids,trees,audited=False):
    for identifier in ids:
        for tree in trees:
            row=dict(fit_input_id=identifier,tree=tree,status=PASS)
            if audited:row['numerical_fit_verified']=True
            yield row


class ScopeTests(unittest.TestCase):
    def setUp(self):
        self.ids={'a','b'};self.trees={'t1','t2','t3','t4','t5'}
        self.fits=list(rows(self.ids,self.trees));self.audits=list(rows(self.ids,self.trees,True))

    def run_scope(self,fits=None,audits=None):
        return full_scope(self.fits if fits is None else fits,self.audits if audits is None else audits,self.trees,self.ids)

    def test_full_production_sized_cartesian_grid(self):
        ids={str(n) for n in range(75070)}
        production,flags,errors=full_scope(rows(ids,self.trees),rows(ids,self.trees,True),self.trees,ids)
        self.assertEqual(len(production),375350);self.assertEqual(flags,[]);self.assertEqual(errors,[])

    def test_flag_and_error_separately_retained(self):
        self.fits[0]['status']=self.audits[0]['status']=FLAG
        self.fits[1]['status']=self.audits[1]['status']=ERROR;self.audits[1]['numerical_fit_verified']=False
        _,flags,errors=self.run_scope();self.assertEqual(len(flags),1);self.assertEqual(len(errors),1)

    def test_missing_fit_rejected(self):
        with self.assertRaises(ValueError):self.run_scope(fits=self.fits[:-1])

    def test_duplicate_fit_rejected(self):
        with self.assertRaises(ValueError):self.run_scope(fits=self.fits+[self.fits[0]])

    def test_duplicate_audit_rejected(self):
        with self.assertRaises(ValueError):self.run_scope(audits=self.audits+[self.audits[0]])

    def test_unknown_tree_rejected(self):
        self.fits[0]['tree']='wrong'
        with self.assertRaises(ValueError):self.run_scope()

    def test_status_mismatch_rejected(self):
        self.audits[0]['status']=FLAG
        with self.assertRaises(ValueError):self.run_scope()

    def test_unverified_numeric_pass_rejected(self):
        self.audits[0]['numerical_fit_verified']=False
        with self.assertRaises(ValueError):self.run_scope()

    def test_error_not_qualified_as_fit(self):
        self.fits[0]['status']=self.audits[0]['status']=ERROR
        with self.assertRaises(ValueError):self.run_scope()


if __name__=='__main__':unittest.main()
