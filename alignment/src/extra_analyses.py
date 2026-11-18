"""Supplementary analyses cited in the report: popularity and length baselines, neighbor
properties, per-class / per-degree breakdowns, example neighbors, and the auto-ubiquitination
sensitivity check. Writes results/sw_per_query.csv and results/sw_extra_analyses.txt."""
import sys; sys.path.insert(0,"src")
import numpy as np, pandas as pd
from evalcore import *
KS=[1,3,5,10]
idx,pairs,keep=load(); idx["cls"]=idx.apply(coarse_class,axis=1)
sub=pairs.groupby("e3")["sub"].apply(set)
has=[i for i in keep if idx.loc[i,"UniProt accession"] in sub.index]
sets=[sub[idx.loc[i,"UniProt accession"]] for i in has]; S=share_matrix(sets); n=len(has)
deg=np.array([len(s) for s in sets]); cls=idx.loc[has,"cls"].to_numpy(dtype=object)
genes=idx.loc[has,"Gene symbol"].to_numpy(dtype=object); L=idx.loc[has,"Sequence length (aa)"].to_numpy()
out=[]
def p(*a): s=" ".join(str(x) for x in a); print(s); out.append(s)
fmt=lambda d:{k:round(float(v),4) for k,v in d.items()}
# 1. baselines that ignore sequence similarity
sc=np.tile(deg.astype(float),(n,1)); np.fill_diagonal(sc,-np.inf)
pop=np.argsort(-sc,axis=1,kind="stable")[:,:10]
p("POPULARITY baseline (always the most substrate-rich other E3s):",fmt(metrics(pop,S,KS)))
Ld=np.abs(L[:,None]-L[None,:]).astype(float); np.fill_diagonal(Ld,np.inf)
ln=np.argsort(Ld,axis=1,kind="stable")[:,:10]
p("LENGTH-matched baseline (nearest by sequence length):",fmt(metrics(ln,S,KS)))
# 2. what kind of proteins do the neighbors tend to be?
for name,fn in (("sqrt","results/sw_norm_sqrt.npy"),("min","results/sw_norm_min.npy"),("raw","results/sw_raw.npy")):
    sim=np.load(fn)[np.ix_(has,has)]; nb=topk_neighbors(sim,10); m=metrics(nb,S,KS)
    p(f"SW {name}: prec@5={m['prec@5']:.4f} hit@5={m['hit@5']:.4f} | mean #substrates of top-5 neighbors={deg[nb[:,:5]].mean():.1f} (all {deg.mean():.1f}) | mean length of top-5 neighbors={L[nb[:,:5]].mean():.0f} (all {L.mean():.0f})")
# 3. per-query table and breakdowns (primary: sqrt-normalized)
sim=np.load("results/sw_norm_sqrt.npy")[np.ix_(has,has)]; nb=topk_neighbors(sim,10)
sh5=np.take_along_axis(S,nb[:,:5],axis=1)
q=pd.DataFrame({"gene":genes,"cls":cls,"n_substrates":deg,"prec5":sh5.mean(1),"hit5":sh5.any(1),"top1_sim":sim[np.arange(n),nb[:,0]]})
q.to_csv("results/sw_per_query.csv",index=False)
q["deg_bin"]=pd.cut(q.n_substrates,[0,1,2,5,10,1000],labels=["1","2","3-5","6-10",">10"])
q["sim_bin"]=pd.cut(q.top1_sim,[0,0.1,0.2,0.3,0.5,1.01])
def table(col):
    g=q.groupby(col,observed=True).agg(prec5=("prec5","mean"),hit5=("hit5","mean"),n=("prec5","size")).round(3); return g
p("\nBy number of substrates of the query:\n"+table("deg_bin").to_string())
p("\nBy class (n shown; classes with n<10 are noisy):\n"+table("cls").sort_values("n",ascending=False).to_string())
p("\nBy similarity of the query's closest neighbor:\n"+table("sim_bin").to_string())
# 4. example neighbors
p("")
for g in ["MDM2","FBXW7","BTRC","SMURF1","TRIM25","PRKN","CBL"]:
    w=np.where(genes==g)[0]
    if len(w)==0: continue
    i=w[0]
    desc=", ".join(f"{genes[j]}({sim[i,j]:.2f}, {len(sets[i]&sets[j])} shared)" for j in nb[i,:5])
    p(f"{g} [{cls[i]}, {deg[i]} substrates]: {desc}")
# 5. sensitivity: drop auto-ubiquitination (E3 listed as its own substrate)
p("")
for label,pp in (("with auto-ubiquitination pairs (main)",pairs),("auto-ubiquitination (E3==substrate) removed",pairs[pairs.e3!=pairs["sub"]])):
    sb=pp.groupby("e3")["sub"].apply(set)
    h=[i for i in keep if idx.loc[i,"UniProt accession"] in sb.index]
    st=[sb[idx.loc[i,"UniProt accession"]] for i in h]; S2=share_matrix(st); d2=np.array([len(s) for s in st]); m2=len(h)
    bins=np.digitize(d2,[2,3,6,11,31]); s2=np.load("results/sw_norm_sqrt.npy")[np.ix_(h,h)]; nb2=topk_neighbors(s2,10)
    o=metrics(nb2,S2,[5])["prec@5"]; rng=np.random.default_rng(0); dist=[]
    for _ in range(1000):
        pm=np.arange(m2)
        for b in np.unique(bins):
            w=np.where(bins==b)[0]; pm[w]=rng.permutation(w)
        dist.append(S2[pm[:,None],pm[nb2[:,:5]]].mean())
    dist=np.array(dist); p(f"{label}: E3s={m2}, pairs={len(pp)}, prec@5={o:.4f}, popularity-controlled null={dist.mean():.4f}, z={(o-dist.mean())/dist.std():.1f}")
p(f"(auto-ubiquitination pairs in the evaluation set: {(pairs.e3==pairs['sub']).sum()})")
open("results/sw_extra_analyses.txt","w").write("\n".join(out)+"\n")