"""All-vs-all Smith-Waterman (local) alignment of human E3-ome reference sequences.
Settings (fixed before any evaluation): BLOSUM62, gap open 12 / extend 1 (parasail convention:
gap of length k costs 12 + (k-1)*1, i.e. BLAST's 11 + k), SW local, 32-bit scores (no saturation).
Primary normalization: score / sqrt(self_i * self_j). Sensitivity: score / min(self_i, self_j)."""
import numpy as np, pandas as pd, parasail, time
from multiprocessing import Pool
r=pd.read_excel("data/E3_substrate_pairs_and_e3_roster.xlsx",sheet_name="E3 roster").dropna(subset=["Protein sequence"]).reset_index(drop=True)
seqs=r["Protein sequence"].tolist(); n=len(seqs)
OPEN,EXT=12,1
def row(i):
    return i,[parasail.sw_striped_32(seqs[i],seqs[j],OPEN,EXT,parasail.blosum62).score for j in range(i,n)]
if __name__=="__main__":
    t=time.time()
    raw=np.zeros((n,n),dtype=np.float64)
    with Pool(2) as p:
        for i,vals in p.imap_unordered(row,range(n),chunksize=8):
            raw[i,i:]=vals
    raw=np.triu(raw)+np.triu(raw,1).T
    self_=np.diag(raw).copy()
    sq=raw/np.sqrt(np.outer(self_,self_)); mn=raw/np.minimum.outer(self_,self_)
    np.save("results/sw_raw.npy",raw); np.save("results/sw_norm_sqrt.npy",sq); np.save("results/sw_norm_min.npy",mn)
    r[["Gene symbol","UniProt accession","E3-ome class","E3-ome role","Receptor family","Sequence length (aa)"]].assign(self_score=self_).to_csv("results/e3_index.csv",index=False)
    print(f"{n} proteins, {n*(n-1)//2} pairs, {time.time()-t:.0f}s")
    print("self-score range",self_.min(),self_.max(),"| norm_sqrt off-diag range",np.round([sq[np.triu_indices(n,1)].min(),sq[np.triu_indices(n,1)].max()],3))
