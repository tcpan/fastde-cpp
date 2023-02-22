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
#transpose it.
mt = m.transpose().tocsc()

h5_io.write_h5_spmat("../data/test_spmat_csc_transposed_csc.h5", mt)

# %%
#transpose to dense
mtd = mt.todense().reshape(mt.shape, order='F')

h5_io.write_h5_matrix("../data/test_spmat_csc_transpose_to_dense_c.h5", mtd)

# %%



# %%
# create transpose from spmat
m = h5_io.read_h5_matrix("../data/test_matrix_c.h5")
# print(repr(m))

# %%
#transpose to dense
mr = m.reshape(m.shape, order='C')

h5_io.write_h5_matrix("../data/temp_matrix_c_r.h5", mr)

# %%
m = h5_io.read_h5_matrix("../data/test_matrix_r.h5")
# print(repr(m))

# %%
#transpose to dense
mc = m.reshape(m.shape, order='F')

h5_io.write_h5_matrix("../data/temp_matrix_r_c.h5", mc)

h5_io.write_h5_matrix("../data/temp_matrix_r.h5", m)


# %%
mt = m.transpose()
h5_io.write_h5_matrix("../data/temp_matrix_r_transpose.h5", mt)
