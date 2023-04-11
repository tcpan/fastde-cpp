#%%
# loading the h5_io package.
import sys
sys.path.append(".")
import h5_io

from scipy import sparse
import numpy as np

import h5py


# %%
# import os
# data_dir = os.path.realpath(os.path.dirname(os.path.realpath(__file__)) + "/../data")


# %%
# create transpose from spmat
m = h5_io.read_h5_spmat("../data/test_spmat_csc.h5", "csc")
scale_factor = 0.5
rows = m.shape[0]
cols = m.shape[1]


# %%
# create a label vector
num_labels = 30
labels = np.random.randint(1, num_labels+1, rows)

h5_io.write_h5_vector("../data/test_label_vec.h5", labels)


# %%
# compute wilcoxon.
import scipy.stats as stats

# for each label, compute in vs out classes.

# gene by gene.  not sure if mannwhitneyu can handel array.
pv = np.zeros(shape = [num_labels, cols], dtype = np.float64, order='F')
pv2 = np.zeros(shape = [num_labels, cols], dtype = np.float64, order='C')
for g in range(cols):
    gene_data = m[:, g].todense()  # sparse matrix, so no squeeze.  
    for l in range(1, num_labels+1):
        sample1 = gene_data[np.where(labels == l)[0].tolist()]
        sample2 = gene_data[np.where(labels != l)[0].tolist()]
        # isnan is not supported by sparse matrix in stats.  so need to convert to dense.
        _, pval = stats.mannwhitneyu(sample1, sample2, use_continuity=True, alternative='two-sided', method='asymptotic')
        pv[l-1, g] = pval[0, 0]
        pv2[l-1, g] = pval[0, 0]

h5_io.write_h5_matrix("../data/test_spmat_wilcox_pval_matrix_c.h5", pv)
h5_io.write_h5_matrix("../data/test_spmat_wilcox_pval_matrix_r.h5", pv2)



# %%
