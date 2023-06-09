#%%
# loading the h5_io package.
import sys
sys.path.append(".")
import h5_io

from scipy import sparse
import numpy as np

import math

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
if (os.path.isfile("../data/test_label_vec.h5")):
    labels = h5_io.read_h5_vector("../data/test_label_vec.h5")

    unique_labs = set(labels)
    num_labels = len(unique_labs)

else:
    # create a label vector
    num_labels = 30
    labels = np.random.randint(1, num_labels+1, rows)

    h5_io.write_h5_vector("../data/test_label_vec.h5", labels)


# %%
# compute fc.
import scipy.stats as stats

# for each label, compute in vs out classes.

# gene by gene.  not sure if mannwhitneyu can handel array.
fc = np.zeros(shape = [num_labels, cols], dtype = np.float64, order='F')
p1 = np.zeros(shape = [num_labels, cols], dtype = np.float64, order='F')
p2 = np.zeros(shape = [num_labels, cols], dtype = np.float64, order='F')
for g in range(cols):
    gene_data = m[:, g].todense()  # sparse matrix, so no squeeze.  
    for l in range(1, num_labels+1):
        sample1 = gene_data[np.where(labels == l)[0].tolist()]
        sample2 = gene_data[np.where(labels != l)[0].tolist()]
        # isnan is not supported by sparse matrix in stats.  so need to convert to dense.
        thresh_min = 0.25
        fc[l-1, g] = math.log2(sum(sample1)/len(sample1) + 1) - math.log2(sum(sample2)/len(sample2) + 1)
        p1[l-1, g] = sum(sample1 > thresh_min) / len(sample1)
        p2[l-1, g] = sum(sample2 > thresh_min) / len(sample2)

h5_io.write_h5_matrix("../data/test_spmat_foldchange_matrix_c.h5", fc)
h5_io.write_h5_matrix("../data/test_spmat_foldchange_p1_matrix_c.h5", p1)
h5_io.write_h5_matrix("../data/test_spmat_foldchange_p2_matrix_c.h5", p2)



# %%
