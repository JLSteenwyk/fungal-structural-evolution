"""Orientation-consistent metadata differences; no structural outcome enters scoring."""
import math
BANDS={'tight':(1.1,5.,.05),'moderate':(1.25,10.,.1),'wide':(1.5,15.,.2)}

def differences(target,background):
    out={}
    for order,sides in enumerate([('a','b'),('b','a')]):
        lr=[];pc=[];lc=[]
        for ts,bs in zip(['a','b'],sides):
            tl=float(target['length_'+ts]);bl=float(background['length_'+bs])
            tp=float(target['mean_ca_plddt_'+ts]);bp=float(background['mean_ca_plddt_'+bs])
            tf=float(target['fraction_ca_plddt_below50_'+ts]);bf=float(background['fraction_ca_plddt_below50_'+bs])
            if not all(math.isfinite(x) for x in [tl,bl,tp,bp,tf,bf]) or min(tl,bl)<=0 or not all(0<=x<=100 for x in [tp,bp]) or not all(0<=x<=1 for x in [tf,bf]):raise ValueError('Invalid endpoint covariates')
            lr.append(max(tl/bl,bl/tl));pc.append(abs(tp-bp));lc.append(abs(tf-bf))
        out['length_ratio_'+str(order)]=max(lr);out['plddt_difference_'+str(order)]=max(pc);out['lowconf_difference_'+str(order)]=max(lc)
    for band,(length,confidence,fraction) in BANDS.items():
        out['within_'+band]=int(any(out['length_ratio_'+str(o)]<=length and out['plddt_difference_'+str(o)]<=confidence and out['lowconf_difference_'+str(o)]<=fraction for o in [0,1]))
    out['shared_genes']=len({target['gene_a'],target['gene_b']}&{background['gene_a'],background['gene_b']})
    out['shared_models']=len({(target['model_id_'+s],target['version_'+s]) for s in ['a','b']}&{(background['model_id_'+s],background['version_'+s]) for s in ['a','b']})
    out['shared_sequences']=len({target['sequence_sha256_'+s] for s in ['a','b']}&{background['sequence_sha256_'+s] for s in ['a','b']})
    return out
