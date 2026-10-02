"""Reusable full-cohort kernel products; no covariance fitting or basis selection."""
import numpy as np
from scipy import sparse
from scipy.linalg import qr
from covariance_basis_audit import _diagnostics


class ComponentKernelProducts:
    def __init__(self, labels, incidence, diagonal):
        labels=np.asarray(labels); n=len(labels)
        if labels.ndim!=1 or not n:raise ValueError('Nonempty component labels required')
        self.n=n; self.diagonal=np.asarray(diagonal,dtype=float)
        if self.diagonal.shape!=(n,) or not np.isfinite(self.diagonal).all() or np.any(self.diagonal<=0):
            raise ValueError('Positive finite residual diagonal required')
        if any(not isinstance(name,str) or not name for name in incidence) or {'residual','species'}&set(incidence):
            raise ValueError('Nonempty entity names must avoid reserved names')
        _,codes=np.unique(labels,return_inverse=True)
        order=np.argsort(codes,kind='stable'); self.parts=np.split(order,np.flatnonzero(np.diff(codes[order]))+1)
        self.operators={}; self.local=[]
        for name,source in incidence.items():
            z=sparse.csr_matrix(source,dtype=float,copy=True)
            if z.shape[0]!=n or not np.isfinite(z.data).all():raise ValueError('Invalid entity incidence')
            z.sum_duplicates(); z.eliminate_zeros(); z.sort_indices()
            # Remove unused columns once; mathematical loadings remain unchanged.
            z=z[:,np.unique(z.indices)].tocsr()
            rows,columns=z.nonzero(); lo=np.full(z.shape[1],n); hi=np.full(z.shape[1],-1)
            np.minimum.at(lo,columns,codes[rows]); np.maximum.at(hi,columns,codes[rows])
            if np.any(lo!=hi):raise ValueError('Entity crosses declared components')
            self.operators[name]=z
        self.names=['residual',*self.operators,'species']; self.raw=np.zeros((len(self.names),len(self.names)))
        self.raw[0,0]=self.diagonal@self.diagonal
        for part in self.parts:
            local=[z[part][:,np.unique(z[part].indices)].tocsr() for z in self.operators.values()]
            self.local.append(local)
            kernels=[(z@z.T).toarray() for z in local]
            if kernels:
                flat=np.asarray(kernels).reshape(len(kernels),-1)
                self.raw[1:-1,1:-1]+=flat@flat.T
                for i,kernel in enumerate(kernels,1):self.raw[0,i]+=self.diagonal[part]@np.diag(kernel)
        self.raw[1:,0]=self.raw[0,1:]

    def tree(self, factor):
        return TreeKernelProducts(self,factor)


class TreeKernelProducts:
    def __init__(self,bank,factor):
        self.bank=bank; self.factor=np.asarray(factor,dtype=float)
        if self.factor.ndim!=2 or self.factor.shape[0]!=bank.n or not np.isfinite(self.factor).all():
            raise ValueError('Finite n-row species factor required')
        self.raw=bank.raw.copy(); self.raw[-1,-1]=np.sum((self.factor.T@self.factor)**2)
        self.raw[0,-1]=self.raw[-1,0]=bank.diagonal@np.sum(self.factor*self.factor,axis=1)
        for part,local in zip(bank.parts,bank.local):
            for i,z in enumerate(local,1):
                cross=z.T@self.factor[part]; self.raw[i,-1]+=np.sum(cross*cross)
        self.raw[-1,1:-1]=self.raw[1:-1,-1]

    def audit(self, design):
        x=np.asarray(design,dtype=float); bank=self.bank
        if x.ndim!=2 or x.shape[0]!=bank.n or not 0<x.shape[1]<bank.n or not np.isfinite(x).all():
            raise ValueError('Finite active design with residual dimensions required')
        maximum=np.max(abs(x),axis=0)
        if np.any(maximum==0):raise ValueError('Inactive zero design columns require separate accounting')
        scaled=x/maximum; scaled/=np.sqrt(np.sum(scaled*scaled,axis=0))
        singular=np.linalg.svd(scaled,compute_uv=False)
        if singular[-1]<=10*max(x.shape)*np.finfo(float).eps*singular[0]:
            raise ValueError('Fixed design rank or boundary requires review')
        q,_=qr(scaled,mode='economic'); condition=float(singular[0]/singular[-1])
        images=[bank.diagonal[:,None]*q]
        images.extend(z@(z.T@q) for z in bank.operators.values())
        images.append(self.factor@(self.factor.T@q)); images=np.asarray(images)
        flat=images.reshape(len(images),-1); pair=flat@flat.T
        cores=np.asarray([q.T@value for value in images]).reshape(len(images),-1); small=cores@cores.T
        projected=self.raw-2*pair+small
        dimension=max(bank.n,self.factor.shape[1],x.shape[1],*(z.shape[1] for z in bank.operators.values()),1)
        envelope=64*np.finfo(float).eps*dimension
        raw_error=envelope*abs(self.raw)
        projected_error=envelope*condition*(abs(self.raw)+2*abs(pair)+abs(small))
        return dict(kernel_names=bank.names,records=bank.n,fixed_effect_columns=x.shape[1],
            residual_dimension=bank.n-x.shape[1],family_components=len(bank.parts),
            normalized_design_condition_number=condition,raw_gram=self.raw.tolist(),projected_gram=projected.tolist(),
            raw_roundoff_envelope=raw_error.tolist(),projected_roundoff_envelope=projected_error.tolist(),
            exactly_zero_incidence_names=[name for name,z in bank.operators.items() if not z.nnz],
            raw_diagnostics=_diagnostics(self.raw,raw_error,bank.names),
            reml_diagnostics=_diagnostics(projected,projected_error,bank.names),covariance_model_selected=False,
            scope='Reusable component and species kernel products, with full fixed-design residual-space '
                  'qualification. Empty entity columns are removed without changing kernels. '
                  'No variance fit, automatic basis selection, clipping or inferential calibration.')
