from scipy.stats import entropy
import numpy as np
import operator
from scipy.spatial import distance
import pandas as pd
import warnings
import gseapy as gp
from scipy import sparse
import scanpy as sc

warnings.filterwarnings('ignore')

def score_signatures(adata, markers):
    for sig, genes in markers.items():
        sc.tl.score_genes(adata,gene_list=sig,score_name=sig)

def get_expression(adata, gene):
    return adata.X[:,adata.var.index.tolist().index(gene)].T.todense().tolist()[0]

def deg_df(adata, category, adjp_thresh=0.001, lfc_thresh=0.25):
    sc.tl.rank_genes_groups(adata,category)
    categories = set(adata.obs[category])
    dfs = []
    for cat in categories:
        df = sc.get.rank_genes_groups_df(adata,cat)
        df[category] = cat 
        dfs.append(df)
    df = pd.concat(dfs)
    df = df[df["pvals_adj"] < adjp_thresh]
    df = df[df["logfoldchanges"].abs() > lfc_thresh]
    return df

def gene_linear_regression(adata, key, significance=0.001):
    cts = []
    pvs = []
    genes = []
    coeffs = []
    df = adata.obs.copy()
    key_values = set(adata.ob[key].tolist())
    for gene in tqdm.tqdm(adata.var.index.tolist()):
        df["Expression"] = adata.X[:,adata.var.index.tolist().index(gene)].T.todense().tolist()[0]
        model = smf.ols('Expression ~ {}'.format(key), data=df).fit()
        coef = model.params['Timepoint[T.{}]'.format(key_values[-1])]  # Coefficient for the T2 category compared to T1 (base category)
        p_value = model.pvalues['Timepoint[T.{}]'.format(key_values[-1])]  # P-value for the T2 category compared to T1 (base category)
        coeffs.append(coef)
        pvs.append(p_value)
        genes.append(gene)
    df = pandas.DataFrame.from_dict({"Gene":genes,"Coefficient":coeffs, "P-value":pvs})
    return df[df["P-value"] < significance]