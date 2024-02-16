from scipy.stats import entropy
import statsmodels.formula.api as smf
import statsmodels.api as sm
import numpy as np
import pandas as pd
import warnings
import gseapy as gp
from scipy import sparse
import scanpy as sc
import collections
import tqdm
import scipy

warnings.filterwarnings('ignore')


def normalized_exponential_vector(values, temperature=0.000001):
    assert temperature > 0, "Temperature must be positive"
    exps = np.exp(values / temperature)
    return exps / np.sum(exps)

def classify(adata, markers, temperature=0.001):
    scores = []
    for ph, genes in markers.items():
        sc.tl.score_genes(adata,score_name="{}_SCORE".format(ph),gene_list=genes)
        scores.append("{}_SCORE".format(ph))
    mat = adata.obs[scores].to_numpy()
    cts = []
    probabs = collections.defaultdict(list)
    found_nan = False
    for x in mat:
        probs = normalized_exponential_vector(x, temperature=temperature)
        if np.isnan(probs).any():
            found_nan = True
            uniform_p = 1. / len(probs)
            for ph, p in zip(scores,probs):
                probabs[ph.replace("_SCORE"," Pseudo-probability")].append(uniform_p)
            ct = scores[0].replace("_SCORE","")
        else:
            for ph, p in zip(scores,probs):
                probabs[ph.replace("_SCORE"," Pseudo-probability")].append(p)
            ct = scores[np.argmax(probs)].replace("_SCORE","")
        cts.append(ct)
    if found_nan:
        message = "Some cells had uniform phenotype probabilities, resulting in a random phenotype classification. You may want to increase the temperature."
        print(message)
        warnings.warn(message, RuntimeWarning)
    for n, p in probabs.items():
        adata.obs[n] = p
    adata.obs['genevector'] = cts
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

def gene_linear_regression(adata, key, coeff, significance=0.001):
    cts = []
    pvs = []
    genes = []
    coeffs = []
    df = adata.obs.copy()
    key_values = set(adata.ob[key].tolist())
    key_coeff = key_values.index(coeff)
    for gene in tqdm.tqdm(adata.var.index.tolist()):
        df["Expression"] = adata.X[:,adata.var.index.tolist().index(gene)].T.todense().tolist()[0]
        model = smf.ols('Expression ~ {}'.format(key), data=df).fit()
        coef = model.params['Timepoint[T.{}]'.format(key_coeff)]  
        p_value = model.pvalues['Timepoint[T.{}]'.format(key_coeff)]
        coeffs.append(coef)
        pvs.append(p_value)
        genes.append(gene)
    df = pd.DataFrame.from_dict({"Gene":genes,"Coefficient":coeffs, "P-value":pvs})
    return df[df["P-value"] < significance]

# def correlation(adata)
#     mis = []
#     pids  = []
#     conds = []
#     correlation = []
#     gt = []
#     for condition in set(adata.obs['subtype']):
#         fdata = adata[adata.obs["subtype"]==condition]
#         mat = []
#         index = []
#         for ct in tqdm.tqdm(set(fdata.obs['pid'])):
#             xdata = fdata[fdata.obs["pid"]==ct]
#             batf = xdata.X[:,xdata.var.index.tolist().index("KLF2")]
#             ncells = len(xdata.obs.index)
#             for gene in ["LEF1","SELL","CCR7","TCF7","IL7R","GZMK","LAG3","GZMA","GZMB"]:
#                 target = xdata.X[:,xdata.var.index.tolist().index(gene.upper())]
#                 if len(target) == 1: continue
#                 x = []
#                 y = []
#                 for c,z in zip(batf, target):
#                     if c > 0 or z > 0:
#                         x.append(c)
#                         y.append(z)
#                 pxy, xedges, yedges = numpy.histogram2d(x, y, density=True)
#                 pxy = pxy / pxy.sum()
#                 px = np.sum(pxy, axis=1)
#                 px = px / px.sum()
#                 py = np.sum(pxy, axis=0)
#                 py = py / py.sum()
#                 px_py = px[:, None] * py[None, :]
#                 nzs = pxy > 0
#                 mi = np.sum(pxy[nzs] * numpy.log2((pxy[nzs] / px_py[nzs])))
#                 c = scipy.stats.pearsonr(target,batf)
#                 mis.append(mi)
#                 correlation.append(c.statistic)
#                 conds.append(condition)
#                 pids.append(ct)
#                 gt.append(gene)
#     dfx = pd.DataFrame.from_dict({"Patient":pids,"Condition":conds,"r2":correlation,"MI":mis,"Target":gt})
#     dfx