import sys; sys.path.insert(0,"src")
import numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from evalcore import *
BLUE,ORANGE,GRAY,INK,INK2,GRID="#2a78d6","#eb6834","#8c8b86","#0b0b0b","#52514e","#e6e5e1"
plt.rcParams.update({"font.family":"DejaVu Sans","font.size":10,"axes.edgecolor":GRID,"axes.labelcolor":INK2,"xtick.color":INK2,"ytick.color":INK2,"axes.spines.top":False,"axes.spines.right":False,"axes.grid":True,"grid.color":GRID,"grid.linewidth":0.8,"axes.axisbelow":True})
idx,pairs,keep=load(); sub=pairs.groupby("e3")["sub"].apply(set)
has=[i for i in keep if idx.loc[i,"UniProt accession"] in sub.index]
sets=[sub[idx.loc[i,"UniProt accession"]] for i in has]; S=share_matrix(sets); n=len(has)
deg=np.array([len(s) for s in sets]); bins=np.digitize(deg,[2,3,6,11,31])
sim=np.load("results/sw_norm_sqrt.npy")[np.ix_(has,has)]; nb=topk_neighbors(sim,10)
KS=[1,3,5,10]; obs=[np.take_along_axis(S,nb[:,:k],axis=1).mean() for k in KS]
sc=np.tile(deg.astype(float),(n,1)); np.fill_diagonal(sc,-np.inf); pop=np.argsort(-sc,axis=1,kind="stable")[:,:10]
popv=[np.take_along_axis(S,pop[:,:k],axis=1).mean() for k in KS]
rng=np.random.default_rng(0); G=[];D=[]
for _ in range(1000):
    p=rng.permutation(n); G.append([S[p[:,None],p[nb[:,:k]]].mean() for k in KS])
    p=np.arange(n)
    for b in np.unique(bins):
        m=np.where(bins==b)[0]; p[m]=rng.permutation(m)
    D.append([S[p[:,None],p[nb[:,:k]]].mean() for k in KS])
G=np.array(G);D=np.array(D)
fig,ax=plt.subplots(1,2,figsize=(11,4.4),gridspec_kw={"width_ratios":[1,1.15]})
a=ax[0]; x=np.arange(len(KS))
a.fill_between(x,np.percentile(D,2.5,0),np.percentile(D,97.5,0),color=GRAY,alpha=.28,lw=0,label="Popularity-controlled shuffle (95% band)")
a.plot(x,D.mean(0),color=GRAY,lw=1.5,ls=":")
a.plot(x,G.mean(0),color=GRAY,lw=1.2,ls="--",label="Simple label shuffle (mean)")
a.plot(x,popv,color=ORANGE,lw=2,marker="o",ms=7,mec="#fcfcfb",mew=2,label="Always return most-studied E3s")
a.plot(x,obs,color=BLUE,lw=2,marker="o",ms=7,mec="#fcfcfb",mew=2,label="Smith-Waterman nearest neighbors")
a.set_xticks(x); a.set_xticklabels([f"k={k}" for k in KS]); a.set_ylabel("Precision@k\n(top-k neighbors sharing a substrate)")
a.set_ylim(0,0.44); a.set_title("a  Neighbors vs. baselines",loc="left",color=INK,fontsize=11,fontweight="bold")
a.text(2.08,obs[2]-0.018,f"{obs[2]:.2f}",color=INK,fontsize=9,va="top"); a.text(2.08,popv[2]+0.008,f"{popv[2]:.2f}",color=INK,fontsize=9)
a.legend(loc="upper left",frameon=False,fontsize=8.2,ncol=1)
b=ax[1]; d=pd.read_csv("results/sw_by_degree_bin.csv"); y=np.arange(len(d))
w=0.34
b.bar(y-w/2,d.prec5_null,w-0.04,color=GRAY,alpha=.55,label="Popularity-controlled shuffle")
b.errorbar(y-w/2,d.prec5_null,yerr=[d.prec5_null-d.null_lo,d.null_hi-d.prec5_null],fmt="none",ecolor=INK2,lw=1,capsize=3)
b.bar(y+w/2,d.prec5_obs,w-0.04,color=BLUE,label="Smith-Waterman neighbors")
for i,r in d.iterrows(): b.text(i+w/2,r.prec5_obs+0.012,f"{r.prec5_obs:.2f}",ha="center",fontsize=8.5,color=INK)
b.set_xticks(y); b.set_xticklabels([f"{g}\n(n={n_})" for g,n_ in zip(d.degree_bin,d.n)]); b.set_xlabel("Substrates documented for the query E3")
b.set_ylabel("Precision@5"); b.set_ylim(0,0.68); b.grid(axis="x",visible=False)
b.set_title("b  Observed vs. null, by how well-studied the query is",loc="left",color=INK,fontsize=11,fontweight="bold")
b.legend(loc="upper left",frameon=False,fontsize=8.5)
fig.tight_layout(); fig.savefig("results/fig1_sw_vs_baselines.png",dpi=200,facecolor="#fcfcfb")
# Fig 2: sensitivity (prec@5 observed vs matching global-shuffle null)
e=pd.read_csv("results/sw_evaluation.csv")
lab={"MAIN: SW sqrt-norm, all scoreable E3s":"Main (all 374 E3s)","Within-class shuffle (labels permuted inside class)":"Main, within-class shuffle null","Normalization variant: min":"Normalize by min self-score","Normalization variant: raw":"Raw score (not normalized)","Set aside RBX1, RNF7, ANAPC11":"Set aside RBX1, RNF7, ANAPC11","Set aside 10 most-studied E3s (by #substrates)":"Set aside 10 most-studied E3s","Exclude close paralogs (norm SW >= 0.5)":"Exclude paralogs, score ≥ 0.5 (37 pairs)","Exclude close paralogs (norm SW >= 0.3)":"Exclude paralogs, score ≥ 0.3 (121 pairs)","Exclude close paralogs (norm SW >= 0.2)":"Exclude paralogs, score ≥ 0.2 (286 pairs)"}
e["lab"]=e.setting.map(lab); e=e.iloc[::-1].reset_index(drop=True)
f,c=plt.subplots(figsize=(8.6,4.6)); yy=np.arange(len(e))
c.hlines(yy,e["prec@5_null_lo"],e["prec@5_null_hi"],color=GRAY,lw=6,alpha=.45,label="Simple label shuffle (95% band)")
c.plot(e["prec@5"],yy,"o",color=BLUE,ms=8,mec="#fcfcfb",mew=2,label="Observed precision@5",zorder=3)
for i,r in e.iterrows(): c.text(r["prec@5"]+0.006,i,f"{r['prec@5']:.3f}",va="center",fontsize=8.5,color=INK)
c.set_yticks(yy); c.set_yticklabels(e.lab); c.set_xlabel("Precision@5 (Smith-Waterman, sqrt-normalized unless stated)"); c.set_xlim(0,0.2); c.grid(axis="y",visible=False)
c.legend(loc="lower center",bbox_to_anchor=(0.4,1.0),ncol=2,frameon=False,fontsize=8.5); f.tight_layout(); f.savefig("results/fig2_sensitivity.png",dpi=200,facecolor="#fcfcfb")
print("ok",np.round(obs,3),np.round(popv,3))
