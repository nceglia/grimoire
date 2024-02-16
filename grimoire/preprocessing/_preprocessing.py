from scipy.stats import entropy
import scanpy as sc
import scanpy.external as sce
import numpy as np
import tqdm
import scanpy.external as sce
import glob
import os


import warnings
warnings.filterwarnings('ignore')

def load_directory(directory, filetype="h5", key="sample"):
    path = os.path.join(directory,"*.{}".format(filetype))
    adatas = []
    samples = []
    for x in glob.glob(path):
        if filetype == "h5":
            adata = sc.read_10x_h5(x)
        else:
            adata = sc.read(x)
        sample = x.replace(".{}".format(filetype),"")
        adata.var_names_make_unique()
        samples.append(sample)
        adatas.append(adata)
    adata = adatas[0].concatenate(adatas[1:],batch_key=key,batch_categories=samples)
    return adata

def run_pca_workflow(adata):
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    return adata

def remove_meaningless_genes(adata, include_mt=True, include_rp=True, include_mtrn=True, include_hsp=True):
    genes = [x for x in adata.var.index.tolist() if "RIK" not in x.upper()]
    genes = [x for x in genes if "GM" not in x]
    genes = [x for x in genes if "-" not in x or "HLA" in x]
    genes = [x for x in genes if "." not in x or "HLA" in x]
    genes = [x for x in genes if "LINC" not in x.upper()]
    if include_mtrn:
        genes = [x for x in adata.var.index.tolist() if "MTRN" not in x]
    if include_hsp:
        genes = [x for x in adata.var.index.tolist() if "HSP" not in x]
    if include_mt:
        genes = [x for x in genes if "MT-" not in x.upper()]
    if include_rp:
        genes = [x for x in genes if "RP" not in x.upper()]
    adata = adata[:,genes]
    return adata

def run_harmony(adata,batch_key):
    sc.tl.pca(adata)
    sce.pp.harmony_integrate(adata, batch_key)
    sc.pp.neighbors(adata,use_rep="X_pca_harmony")
    sc.tl.umap(adata)
    return adata

def gene_entropy(adata, key_added="entropy"):
    X = adata.X.todense()
    X = np.array(X.T)
    gene_to_row = list(zip(adata.var.index.tolist(), X))
    entropies = []
    for _, exp in tqdm.tqdm(gene_to_row):
        counts = np.unique(exp, return_counts = True)
        entropies.append(entropy(counts[1][1:]))
    adata.var[key_added] = entropies
