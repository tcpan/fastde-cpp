#%%
# loading the h5_io package.
import sys
sys.path.append(".")
import h5_io

from scipy import sparse

import numpy as np
import pandas as pd

import h5py



# %%
# create transpose from spmat
m = h5_io.read_h5_spmat("../data/test_spmat_csc.h5", "csc")
# print(m[0:10, 0:10])
# print(m.shape)
# print(m.dtype)
# print(repr(m))

# %%
#to dense. 
mt = np.zeros(m.shape, dtype = m.dtype, order='F')
m.todense(out = mt)   # column major ordering.
# print(mt[0:10, 0:10])
# print(mt.shape)
# print(mt.flags.c_contiguous)
# print(mt.flags.f_contiguous)
# print(mt.dtype)
# print(repr(mt))



# %%
# write out the transposed matrix.
h5_io.write_h5_matrix("../data/test_spmat_csc_to_dense_c.h5", mt)



# %%
#to dense
mt = m.todense().copy(order='C')   # column major ordering.
# print(mt[0:10, 0:10])
# print(mt.shape)
# print(mt.flags.c_contiguous)
# print(mt.flags.f_contiguous)
# print(mt.dtype)
# print(repr(mt))


# %%
# write out the transposed matrix.
h5_io.write_h5_matrix("../data/test_spmat_csc_to_dense_r.h5", mt)

