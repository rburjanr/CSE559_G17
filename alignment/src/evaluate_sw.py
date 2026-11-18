"""Nearest-neighbor evaluation of Smith-Waterman similarity vs documented substrate sharing."""
import sys; sys.path.insert(0,"src")
import numpy as np, pandas as pd, json
from evalcore import *
KS=[1,3,5,10,20]; NPERM=1000; SEED=0
idx,pairs,keep=load(); idx["cls"]=idx.apply(coarse_class,axis=1)
sub=pairs.groupby("e3")["sub"].apply(set)
has=[i for i in keep if idx.loc[i,"UniProt accession"] in sub.index]
SIMS={"sqrt":np.load("results/sw_norm_sqrt.npy"),"min":np.load("results/sw_norm_min.npy"),"raw":np.load("results/sw_raw.npy")}
SQ=SIMS["sqrt"]
def setup(rows):
    sets=[sub[idx.loc[i,"UniProt accession"]] for i in rows]
    return sets,share_matrix(sets),idx.loc[rows,"cls"].to_numpy(dtype=object)
def perms(n,rng,groups=None):
    out=[]
    for _ in range(NPERM):
        if groups is None: out.append(rng.permutation(n))
        else:
            p=np.arange(n)
            for g in np.unique(groups):
                m=np.where(groups==g)[0]; p[m]=rng.permutation(m)
            out.append(p)
    return out
def null_metrics(nb,S,plist):
    kmax=nb.shape[1]; res={f"prec@{k}":[] for k in KS}; res.update({f"hit@{k}":[] for k in KS})
    for p in plist:
        sh=S[p[:,None],p[nb]]
        for k in KS:
            res[f"prec@{k}"].append(sh[:,:k].mean()); res[f"hit@{k}"].append(sh[:,:k].any(axis=1).mean())
    return {k:np.array(v) for k,v in res.items()}
def run(sim,rows,label,exclude_mask=None,within_class=False,seed=SEED):
    sets,S,cls=setup(rows); s=sim[np.ix_(rows,rows)].copy()
    if exclude_mask is not None: s[exclude_mask[np.ix_(rows,rows)]]=-np.inf
    nb=topk_neighbors(s,max(KS))
    obs=metrics(nb,S,KS)
    rng=np.random.default_rng(seed)
    nm=null_metrics(nb,S,perms(len(rows),rng,cls if within_class else None))
    row={"setting":label,"n_queries":len(rows)}
    for k in KS:
        for m in ("prec","hit"):
            key=f"{m}@{k}"; o=obs[key]; d=nm[key]
            row[key]=o; row[key+"_null"]=d.mean(); row[key+"_null_lo"]=np.percentile(d,2.5); row[key+"_null_hi"]=np.percentile(d,97.5)
            row[key+"_p"]=(1+(d>=o).sum())/(1+NPERM); row[key+"_z"]=(o-d.mean())/d.std()
    # class-only reference: expected precision of a random same-class neighbor
    same=(cls[:,None]==cls[None,:]); np.fill_diagonal(same,False)
    ok=same.sum(1)>0
    row["class_only_prec"]=np.mean([S[i,same[i]].mean() for i in np.where(ok)[0]])
    row["class_only_queries"]=int(ok.sum())
    # fraction of top-k neighbors in same class
    row["same_class@5"]=np.take_along_axis(same,nb[:,:5],axis=1).mean()
    row["same_class@10"]=np.take_along_axis(same,nb[:,:10],axis=1).mean()
    return row,nb
if __name__=="__main__":
    out=[]
    r,nb=run(SQ,has,"MAIN: SW sqrt-norm, all scoreable E3s"); out.append(r)
    np.save("results/sw_main_neighbors.npy",nb)
    r,_=run(SQ,has,"Within-class shuffle (labels permuted inside class)",within_class=True); out.append(r)
    for nm in ("min","raw"):
        r,_=run(SIMS[nm],has,f"Normalization variant: {nm}"); out.append(r)
    # set aside heavily studied connector E3s
    acc={g:idx.index[idx["Gene symbol"]==g][0] for g in ("RBX1","RNF7","ANAPC11")}
    rows=[i for i in has if i not in acc.values()]
    r,_=run(SQ,rows,"Set aside RBX1, RNF7, ANAPC11"); out.append(r)
    # set aside 10 most substrate-rich E3s
    deg={i:len(sub[idx.loc[i,"UniProt accession"]]) for i in has}
    top10=sorted(deg,key=deg.get,reverse=True)[:10]
    rows=[i for i in has if i not in top10]
    r,_=run(SQ,rows,"Set aside 10 most-studied E3s (by #substrates)"); out.append(r)
    # paralog exclusion: drop neighbors with normalized SW >= tau
    for tau in (0.5,0.3,0.2):
        mask=SQ>=tau
        r,_=run(SQ,has,f"Exclude close paralogs (norm SW >= {tau})",exclude_mask=mask); out.append(r)
    df=pd.DataFrame(out); df.to_csv("results/sw_evaluation.csv",index=False)
    pd.set_option("display.width",250); pd.set_option("display.max_columns",30)
    cols=["setting","n_queries","prec@1","prec@5","prec@5_null","prec@5_p","prec@10","prec@10_null","hit@5","hit@5_null","class_only_prec","same_class@5"]
    print(df[cols].round(4).to_string(index=False))
    # how many pairs are 'paralog-like' under each threshold
    iu=np.triu_indices(len(has),1); sub_s=SQ[np.ix_(has,has)][iu]
    for tau in (0.5,0.3,0.2): print(f"scoreable pairs with norm SW >= {tau}: {(sub_s>=tau).sum()}")
