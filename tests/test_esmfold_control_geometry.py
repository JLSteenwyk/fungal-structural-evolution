import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from compare_esmfold_controls import compare_arrays

class ControlGeometryTests(unittest.TestCase):
    def arrays(self):
        x=np.arange(60,dtype=float)
        a=np.column_stack((x,np.sin(x),np.cos(x),np.full(60,90.)))
        return a,np.zeros((60,60))
    def test_rigid_transform_is_zero_change(self):
        a,p=self.arrays();b=a.copy()
        rotation=np.array([[0.,-1.,0.],[1.,0.,0.],[0.,0.,1.]])
        b[:,:3]=b[:,:3]@rotation+np.array([2.,5.,-3.])
        r=compare_arrays(a,b,p,p,70)
        self.assertAlmostEqual(r['ca_superposition_rmsd_angstrom'],0.,places=10)
        self.assertAlmostEqual(r['local_distance_mean_absolute_change_angstrom'],0.,places=10)
        self.assertEqual(r['pae10_local_pairs'],r['local_distance_pairs'])
    def test_reverse_direction_in_either_model_excludes_pair(self):
        a,p=self.arrays();q=p.copy();q[10,0]=20
        baseline=compare_arrays(a,a,p,p,70)
        for left,right in ((q,p),(p,q)):
            r=compare_arrays(a,a,left,right,70)
            self.assertEqual(r['pae10_local_pairs'],baseline['local_distance_pairs']-1)
    def test_low_joint_coverage_remains_missing(self):
        a,p=self.arrays();b=a.copy();b[:11,3]=40
        r=compare_arrays(a,b,p,p,70)
        self.assertEqual(r['status'],'insufficient_coverage')
        self.assertEqual(r['matched_residues'],49)
        self.assertEqual(r['ca_superposition_rmsd_angstrom'],'')
    def test_invalid_confidence_rejected(self):
        a,p=self.arrays();a[0,3]=np.nan
        with self.assertRaises(ValueError):compare_arrays(a,a,p,p,70)

if __name__=='__main__':unittest.main()
