import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from statannot import add_stat_annotation
from gseapy import dotplot
import collections
import operator
import itertools
from adjustText import adjust_text
from ..tools._tools import get_expression
import warnings
warnings.filterwarnings('ignore')

dark2_extended_palette = [
    "#1B9E77",  # Dark2 original color
    "#D95F02",  # Dark2 original color
    "#7570B3",  # Dark2 original color
    "#E7298A",  # Dark2 original color
    "#66A61E",  # Dark2 original color
    "#E6AB02",  # Dark2 original color
    "#A6761D",  # Dark2 original color
    "#666666",  # Dark2 original color
    "#2D8E6E",  # Variation
    "#C75505",  # Variation
    "#6C6B9E",  # Variation
    "#D42D7B",  # Variation
    "#5B9730",  # Variation
    "#D4A103",  # Variation
    "#8F6B25",  # Variation
    "#595959",  # Variation
    "#338F65",  # Variation
    "#B74B06",  # Variation
    "#6362C2",  # Variation
    "#C21C6B",  # Variation
    "#4F8E2D",  # Variation
    "#C2A505",  # Variation
    "#865E2F",  # Variation
    "#4D4D4D",  # Variation
    "#3A995C",  # Variation
    "#AC4108",  # Variation
    "#5958D6",  # Variation
    "#B0105C",  # Variation
    "#438529",  # Variation
    "#B0B007"   # Variation
]

sns.set_palette(sns.color_palette(dark2_extended_palette))

def plot_fraction(adata,category,variable,save=None,color=dark2_extended_palette):
    df = adata.obs
    count_df = df.groupby([category, variable]).size().unstack(fill_value=0)
    proportion_df = count_df.divide(count_df.sum(axis=1), axis=0)
    proportion_df.plot(kind="bar", stacked=True, figsize=(10,6),color=color)
    plt.ylabel("Fraction")
    plt.ylim([0, 1])
    plt.legend(title=variable, loc="upper left", bbox_to_anchor=(1, 1))
    plt.tight_layout()
    if save != None:
        plt.savefig(save)

def crosstab_heatmap(adata, variable1, variable2):
    df = adata.obs[[variable1,variable2]]
    df=pd.crosstab(df[variable1],df[variable2],normalize='index')
    return sns.heatmap(df)

def knee_plot(adata, umi_cutoff=1500):
    plt.scatter(
        y = ad.obs.sort_values("n_UMIs", ascending=False)["n_UMIs"],
        x = np.arange(ad.obs.shape[0]),
        s=1, c=ad.obs['percent_mito'])
    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1, 1e6)
    plt.ylim(1, 1e5)
    plt.axhline(y=umi_cutoff)

def boxplot(adata, key, groupby=None, splitby=None, figsize=(5,4),ax=None, order=None, split_order=None, rotation=0, save=None):
    df = adata.obs.copy()
    if key in adata.var.index.tolist():
        exp = get_expression(adata,key)
        df[key] = exp
    if ax == None:
        fig,ax=plt.subplots(1,1,figsize=figsize)
    if order == None:
        order = list(sorted(set(adata.obs[groupby])))
    if splitby != None:
        if split_order != None:
            sns.boxplot(data=df, x=groupby, hue=splitby, y=key,ax=ax, order=order, hue_order=split_order)
        else:
            sns.boxplot(data=df, x=groupby, hue=splitby, y=key,ax=ax, order=order)
    else:
        sns.boxplot(data=df, x=groupby, y=key,ax=ax)
    plt.xticks(rotation=rotation)
    if save != None and ax == None:
        fig.savefig(save)
    return ax

def expression_histogram(adata, gene):
    exp = adata.X[:,adata.var.index.tolist().index(gene)].T.todense().tolist()[0]
    fig,ax = plt.subplots(1,1,figsize=(6,4))
    sns.histplot(exp,ax=ax)
    return ax

def volcano_plot(df,pval_col="pvals_adj",lfc_col='logfoldchanges',gene_col="names",pval_thresh = 0.05,lfc_thresh = 1.0,figsize=(12, 6),top_n=50, title="",max_lfc=10., max_pval=300.):
    df['minus_log10_pvalue'] = -np.log10(df[pval_col])
    df['minus_log10_pvalue'] = np.nan_to_num(df['minus_log10_pvalue'], posinf=max_pval)
    df = df[df[lfc_col].abs() < max_lfc]
    plt.figure(figsize=figsize)
    plt.scatter(df[lfc_col], df['minus_log10_pvalue'], color='gray')
    significant_genes = (df[pval_col] < pval_thresh) & (df[lfc_col].abs() > lfc_thresh)
    plt.scatter(df[lfc_col][significant_genes], df['minus_log10_pvalue'][significant_genes], color='red')
    texts = []
    sorted_df = df[significant_genes].sort_values(by=[pval_col, lfc_col], ascending=[True, False]).head(top_n)
    for i, row in sorted_df.iterrows():
        texts.append(plt.text(row[lfc_col], row['minus_log10_pvalue'], row[gene_col], fontsize=9, ha='right'))
    adjust_text(texts)
    plt.title(title)
    plt.xlabel('Log2(Fold Change)')
    plt.ylabel('-Log10(p-value)')
    plt.axhline(y=-np.log10(pval_thresh), color='black', linestyle='--')
    plt.axvline(x=lfc_thresh, color='black', linestyle='--')
    plt.axvline(x=-lfc_thresh, color='black', linestyle='--')
    plt.tight_layout()
    plt.show()