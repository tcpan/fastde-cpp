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

# nz2 =len(m.nonzero()[1])
# nz1 =len(m.todense().nonzero()[1])   // todense will add entries having the same coordinates.
# print("source matrix has " + str(nz2) + " but " + str(nz1) + " nonzero elements")


# %%
#cbind
mc = m.tocsc()
c_mt = sparse.hstack([mc, mc, mc], format = 'csc')
c_m = c_mt.tocsc()

h5_io.write_h5_spmat("../data/test_spmat_csc_cbind_csc.h5", c_m)


# %%
#rbind - sparse.vstack appears to exclude 0.0 in the input sparse matrix if input format is csc.  fine with csr
# it's likely that sparse.hstack does this for csr.   fine with csc.
# the result is that OBSERVED 0.0 are excluded in these cases.
# for testing, use explicitly convert to the data preserving format. vstack with csr, hstack with csc
# THIS DUE TO X, I, and P value copying and simple transforms.
mr = m.tocsr()
r_mt = sparse.vstack([mr, mr, mr], format = 'csr')
r_m = r_mt.tocsc()
# print(repr(r_m))

h5_io.write_h5_spmat("../data/test_spmat_csc_rbind_csc.h5", r_m)
# %%