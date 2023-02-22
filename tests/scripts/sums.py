#%%
# loading the h5_io package.
import sys
sys.path.append(".")
import h5_io

from scipy import sparse

import h5py



# %%
# create transpose from spmat
m = h5_io.read_h5_spmat("../data/test_spmat_csc.h5", "csc")


# %%
#row sum.  numpy's sum is "along" the specified axis.  along "row" means colsums, not to get sums for each value position of the axis.
rowsum = m.sum(axis = 1)  

h5_io.write_h5_vector("../data/test_spmat_csc_rsum.h5", rowsum)


# %%
#cl sum.   numpy's sum is "along" the specified axis.  along "row" means colsums, not to get sums for each value position of the axis.
colsum = m.sum(axis = 0)

h5_io.write_h5_vector("../data/test_spmat_csc_csum.h5", colsum)
# %%