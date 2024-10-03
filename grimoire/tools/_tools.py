from scipy.stats import entropy
import numpy as np
import pandas as pd
import warnings
import gseapy as gp
from scipy import sparse
import scanpy as sc
import collections

warnings.filterwarnings('ignore')


def normalized_exponential_vector(values, temperature=0.01):
    assert temperature > 0, "Temperature must be positive"
    exps = np.exp(values / temperature)
    return exps / np.sum(exps)

def classify(adata, markers, temperature=0.01, load_probability=False, key="phenotype", probability_suffix="Pseudo-probability"):
    scores = []
    for ph, genes in markers.items():
        sc.tl.score_genes(adata,score_name="{}_SCORE".format(ph),gene_list=genes)
        scores.append("{}_SCORE".format(ph))
    mat = adata.obs[scores].to_numpy()
    cts = []
    probabs = collections.defaultdict(list)
    found_nan = False
    for x in mat:
        probs = normalized_exponential_vector(x,temperature=temperature)
        for ph, p in zip(scores,probs):
            probabs[ph.replace("_SCORE"," {}".format(probability_suffix))].append(p)
        ct = scores[np.argmax(probs)].replace("_SCORE","")
        cts.append(ct)
    adata.uns["probability_columns"] = list(probabs)
    for ph, prob in probabs.items():
        adata.obs[ph] = prob
    adata.obs[key] = cts
    return adata

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

# def gene_linear_regression(adata, key, coeff, significance=0.001):
#     cts = []
#     pvs = []
#     genes = []
#     coeffs = []
#     df = adata.obs.copy()
#     key_values = set(adata.ob[key].tolist())
#     key_coeff = key_values.index(coeff)
#     for gene in tqdm.tqdm(adata.var.index.tolist()):
#         df["Expression"] = adata.X[:,adata.var.index.tolist().index(gene)].T.todense().tolist()[0]
#         model = smf.ols('Expression ~ {}'.format(key), data=df).fit()
#         coef = model.params['Timepoint[T.{}]'.format(key_coeff)]  
#         p_value = model.pvalues['Timepoint[T.{}]'.format(key_coeff)]
#         coeffs.append(coef)
#         pvs.append(p_value)
#         genes.append(gene)
#     df = pd.DataFrame.from_dict({"Gene":genes,"Coefficient":coeffs, "P-value":pvs})
#     return df[df["P-value"] < significance]
