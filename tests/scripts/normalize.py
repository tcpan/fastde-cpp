#%%
# loading the h5_io package.
import sys
sys.path.append(".")
import h5_io

from scipy import sparse
import numpy as np

import h5py



# %%
# create transpose from spmat
m = h5_io.read_h5_spmat("../data/test_spmat_csc.h5", "csc")
scale_factor = 0.5
rows = m.shape[0]


# %%
#compute relative count (el / colsum[k] * scale_factor).  set scale factor to 0.5
colsums = np.squeeze(np.asarray(m.sum(axis = 0)))

colsumMat = np.tile(colsums, (m.shape[0], 1))

out = sparse.csc_matrix(m * scale_factor / colsumMat)

h5_io.write_h5_spmat("../data/test_spmat_csc_relative_count.h5", out)



# %%
#compute lognorm log1p(el / colsum[k] * scale_factor).  set scale factor to 0.5
colsums = np.squeeze(np.asarray(m.sum(axis = 0)))

colsumMat = np.tile(colsums, (rows, 1))

out = np.log1p(sparse.csc_matrix(m * scale_factor / colsumMat))

h5_io.write_h5_spmat("../data/test_spmat_csc_lognorm.h5", out)

# %%
# log1p(x / ( exp( sum(log1p(x[x > 0]), na.rm = TRUE) / length(x) ) ) )
logm = np.log1p(m)
colsums = np.exp(np.squeeze(np.asarray(logm.sum(axis = 0) / rows)))

colsumMat = np.tile(colsums, (rows, 1))

out = np.log1p(sparse.csc_matrix(m / colsumMat))

h5_io.write_h5_spmat("../data/test_spmat_csc_clr.h5", out)
# %%
