import unittest,math
import numpy as np
from Bio.PDB import Structure,Model,Chain,Residue,Atom
from Bio.PDB.SASA import ShrakeRupley

class AccessibilityTests(unittest.TestCase):
    def test_isolated_sphere_and_occlusion(self):
        s=Structure.Structure('test');m=Model.Model(0);s.add(m);c=Chain.Chain('A');m.add(c);r=Residue.Residue((' ',1,' '),'ALA','');c.add(r)
        def atom(name,position,serial):return Atom.Atom(name,np.array(position,dtype=float),90,1,' ',name,serial,element='C')
        r.add(atom('CA',[0,0,0],1));sr=ShrakeRupley(probe_radius=1.4,n_points=960,radii_dict={'C':1.7});sr.compute(m,level='R')
        isolated=4*math.pi*(1.7+1.4)**2;self.assertAlmostEqual(r.sasa,isolated,places=8)
        r.add(atom('CB',[1.5,0,0],2));sr.compute(m,level='R')
        self.assertGreater(r.sasa,isolated);self.assertLess(r.sasa,2*isolated)
        self.assertAlmostEqual(r.sasa,sum(a.sasa for a in r),places=8)

if __name__=='__main__':unittest.main()
