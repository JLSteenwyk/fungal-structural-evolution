"""Independent uniform-residual error-contrast algebra in latent entity space."""
import numpy as np
from scipy import sparse,linalg


class IndependentKernelProducts:
    def __init__(self,labels,incidence):
        self.names=['residual',*incidence,'species']; self.n=len(labels)
        self.parts=[np.flatnonzero(labels==name) for name in sorted(set(labels))]
        self.operators=[z[:,np.unique(z.indices)].tocsr() for z in incidence.values()]
        self.cross={};count=len(incidence)
        self.base=np.zeros((count+2,count+2)); self.base[0,0]=self.n
        self.widths=[len(np.unique(z.indices)) for z in incidence.values()]
        for i,a in enumerate(self.operators):
            self.base[0,i+1]=float(a.multiply(a).sum())
            for j,b in enumerate(self.operators):
                r=(a.T@b).tocsr();self.cross[i,j]=r
                self.base[i+1,j+1]=float(r.multiply(r).sum())
        self.base[1:,0]=self.base[0,1:]

    def tree(self,factor):
        return IndependentTreeProducts(self,factor)


class IndependentTreeProducts:
    def __init__(self,bank,factor):
        self.bank=bank;self.factor=factor;self.raw=bank.base.copy()
        self.species_core=factor.T@factor
        self.raw[0,-1]=self.raw[-1,0]=float(np.sum(factor*factor))
        self.raw[-1,-1]=float(np.sum(self.species_core*self.species_core))
        for i,z in enumerate(bank.operators,1):
            cross=z.T@factor;self.raw[i,-1]=float(np.sum(cross*cross));del cross
        self.raw[-1,1:-1]=self.raw[1:-1,-1]

    def project(self,design):
        bank=self.bank;n=bank.n
        scaled=np.asarray(design,dtype=np.longdouble)
        scaled/=np.max(abs(scaled),axis=0)
        scaled/=np.sqrt(np.sum(scaled*scaled,axis=0))
        scaled=np.asarray(scaled,dtype=float)
        q,_,_=linalg.qr(scaled,mode='economic',pivoting=True)
        singular=linalg.svd(scaled,compute_uv=False,lapack_driver='gesvd')
        assert singular[-1]>10*max(scaled.shape)*np.finfo(float).eps*singular[0]
        size=len(bank.names);p=design.shape[1]
        tiny=np.zeros((size,p,p));correction=np.zeros((size,size))
        species=q.T@self.factor;tiny[-1]=species@species.T
        coefficients=[np.asarray(z.T@q).T for z in bank.operators]
        for i,a in enumerate(coefficients,1):tiny[i]=a@a.T
        for i,a in enumerate(coefficients,1):
            for j,b in enumerate(coefficients,1):
                correction[i,j]=float(np.sum((bank.cross[i-1,j-1].T@a.T).T*b))
            # Associative contraction avoids storing a dense latent-by-species
            # projection for every tree and mode. Dense temporary has n-by-p,
            # while all entity-overlap matrices remain sparse.
            projected_entities=bank.operators[i-1]@a.T
            correction[i,-1]=float(np.sum((projected_entities.T@self.factor)*species))
        correction[-1,1:-1]=correction[1:-1,-1]
        correction[-1,-1]=float(np.sum(self.species_core*(species.T@species)))
        small=np.asarray([[np.sum(a*b) for b in tiny] for a in tiny])
        projected=self.raw-2*correction+small
        projected[0,0]=n-p
        for i in range(1,size):projected[0,i]=projected[i,0]=self.raw[0,i]-np.trace(tiny[i])
        # Supply the same algebraic contraction terms used by the producer's
        # numerical review envelope, derived independently from latent overlaps.
        image_products=correction.copy();image_products[0,0]=p
        core_products=small.copy();core_products[0,0]=p
        for i in range(1,size):
            value=float(np.trace(tiny[i]))
            image_products[0,i]=image_products[i,0]=value
            core_products[0,i]=core_products[i,0]=value
        dimension=max(n,self.factor.shape[1],p,*bank.widths,1)
        envelope=64*np.finfo(float).eps*dimension
        condition=float(singular[0]/singular[-1])
        return dict(raw=self.raw,projected=projected,
            raw_error=envelope*abs(self.raw),
            projected_error=envelope*condition*(abs(self.raw)+2*abs(image_products)+abs(core_products)),
            normalized_design_condition_number=condition)
