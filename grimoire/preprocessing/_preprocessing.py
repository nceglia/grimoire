from scipy.stats import entropy
import numpy as np
import tqdm

import warnings
warnings.filterwarnings('ignore')

def remove_meaningless_genes(adata, include_mt=True, include_rp=True):
    genes = [x for x in adata.var.index.tolist() if "RIK" not in x.upper()]
    genes = [x for x in genes if "GM" not in x]
    genes = [x for x in genes if "-" not in x or "HLA" in x]
    genes = [x for x in genes if "." not in x or "HLA" in x]
    genes = [x for x in genes if "LINC" not in x.upper()]
    if include_mt:
        genes = [x for x in genes if "MT-" not in x.upper()]
    if include_rp:
        genes = [x for x in genes if "RP" not in x.upper()]
    adata = adata[:,genes]
    return adata

def filter_aggresively(adata, min_counts=2, perc_cells=0.1):
    return adata