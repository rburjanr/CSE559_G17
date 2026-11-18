"""Popularity-controlled null: permute substrate-set labels only among E3s in the same substrate-count bin.
This keeps (query degree, neighbor degree) joint structure fixed, so only *which* substrates are shared is randomized."""
import sys; sys.path.insert(0,"src")
import numpy as np, pandas as pd
from evalcore import *
KS=[1,3,5,10]; NPERM=1000
idx,pairs,keep=load(); idx["cls"]=idx.apply(coarse_class,axis=1)
sub=pairs.groupby("e3")["sub"].apply(set)
has=[i for i in keep if idx.loc[i,"UniProt accession"] in sub.index]
sets=[sub[idx.loc[i,"UniProt accession"]] for i in has]; S=share_matrix(sets)
deg=np.array([len(s) for s in sets]); n=len(has)
bins=np.digitize(deg,[2,3,6,11,31])   # 1 | 2 | 3-5 | 6-10 | 11-30 | >30
BL=["1","2","3-5","6-10","11-30",">30"]
def within_bin_perm(rng):
    p=np.arange(n)
    for b in np.unique(bins):
        m=np.where(bins==b)[0]; p[m]=rng.permutation(m)
    return p
res=[]
for name,fn in (("sqrt-norm","results/sw_norm_sqrt.npy"),("min-norm","results/sw_norm_min.npy"),("raw","results/sw_raw.npy")):
    sim=np.load(fn)[np.ix_(has,has)]; nb=topk_neighbors(sim,10); obs=metrics(nb,S,KS)
    rng=np.random.default_rng(0)
    nulls={k:[] for k in obs}; qn={b:[] for b in range(len(BL))}
    shobs=np.take_along_axis(S,nb[:,:5],axis=1).mean(1)
    for _ in range(NPERM):
        p=within_bin_perm(rng); sh=S[p[:,None],p[nb]]
        for k in KS:
            nulls[f"prec@{k}"].append(sh[:,:k].mean()); nulls[f"hit@{k}"].append(sh[:,:k].any(1).mean())
        shq=sh[:,:5].mean(1)
        for b in range(len(BL)): qn[b].append(shq[bins==b].mean())
    for key in ("prec@1","prec@5","hit@5","prec@10"):
        d=np.array(nulls[key]); res.append(dict(sim=name,metric=key,observed=obs[key],null_mean=d.mean(),null_lo=np.percentile(d,2.5),null_hi=np.percentile(d,97.5),z=(obs[key]-d.mean())/d.std(),p=(1+(d>=obs[key]).sum())/(1+NPERM)))
    if name=="sqrt-norm":
        rows=[]
        for b in range(len(BL)):
            d=np.array(qn[b]); o=shobs[bins==b].mean()
            rows.append(dict(degree_bin=BL[b],n=int((bins==b).sum()),prec5_obs=o,prec5_null=d.mean(),null_lo=np.percentile(d,2.5),null_hi=np.percentile(d,97.5),p=(1+(d>=o).sum())/(1+NPERM)))
        pd.DataFrame(rows).to_csv("results/sw_by_degree_bin.csv",index=False)
        print("Per query-degree bin (sqrt-norm), prec@5 observed vs popularity-controlled null:\n",pd.DataFrame(rows).round(3).to_string(index=False))
df=pd.DataFrame(res); df.to_csv("results/sw_degree_controlled_null.csv",index=False)
print("\n",df.round(4).to_string(index=False))
