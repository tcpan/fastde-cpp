
#!python3
#
# Copyright 2020 Georgia Tech Research Corporation
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# Author(s): Tony C. Pan
#

import pandas as pd
import os.path as path
import numpy as np
import h5py
from scipy import sparse

def read_h5_df(filename, dtype=np.float64):
    return pd.read_hdf(filename, 'matrix')

def write_h5_df(df, filename):
    df.to_hdf(filename, 'matrix', mode='w')
    pass


# h5 matix for us should always be row major.
def read_h5_matrix(filename, dtype=np.float64):
    f = h5py.File(filename,'r')
    g = f['matrix']
    # does this always read as row major?
    M = g['block0_values'][:]
    order = g.attrs['order']
    if order == 'F':
        M2 = M.transpose().copy(order = 'F')
    else:
        M2 = M.copy(order ='C')

    f.close()
    return M2

def write_h5_matrix(filename, M):
    # this seems to like to write row by row, rather than via contiguous memory.
    # current fastde code (csc and _c versions) expects column major storage in memory.
    f = h5py.File(filename,'w')
    g = f.create_group('matrix')
    # 
    g.attrs['nblocks'] = 1
    g.attrs['ndim'] = 2
    g.attrs['shape'] = M.shape
    if (M.flags.f_contiguous):
        g.attrs['order'] = 'F'
        # convert to transpose to write out as column major.
        Mt = M.transpose()
        g.create_dataset('block0_values', data=Mt)
    else:
        g.attrs['order'] = 'C'

        # it seems c_contiguous or f_contiguous have no effect.  output is always treated as row major.
        g.create_dataset('block0_values', data=M)

    f.close()
    pass

# input is a scipy spmat
def write_h5_spmat(filename, M):
    f = h5py.File(filename,'w')
    g = f.create_group('sparse_matrix')
    g.create_dataset('x',data=M.data.squeeze())
    g.create_dataset('p',data=M.indptr.squeeze())
    g.create_dataset('i',data=M.indices.squeeze())
    g.attrs['nblocks'] = 1
    g.attrs['ndim'] = 2
    g.attrs['format'] = M.getformat()
    g.attrs['shape'] = M.shape
    f.close()
    pass

# allowed format:  csc, csr
def read_h5_spmat(filename, format = 'csc'):
    f = h5py.File(filename,'r')
    print(list(f.keys()))
    print(list(f['sparse_matrix'].keys()))

    g = f['sparse_matrix']
    M = sparse.csc_matrix((g['x'][:],g['i'][:], g['p'][:]), g.attrs['shape'])

    if format == g.attrs['format']:
        M2 = M
    else:
        if format == 'csc':
            M2 = M.tocsc()
        else:
            M2 = M.tocsr()

    f.close()

    return M2

def write_h5_vector(filename, V):
    Va = np.array(V).squeeze()  # input may be a np matrix (e.g. from sum() function)
    f = h5py.File(filename,'w')
    g = f.create_group('vector')
    g.create_dataset('block0_values',data=Va)
    g.attrs['nblocks'] = 1
    g.attrs['ndim'] = 1
    g.attrs['shape'] = Va.shape
    f.close()
    pass

def read_h5_vector(filename):
    f = h5py.File(filename,'r')
    print(list(f.keys()))
    print(list(f['vector'].keys()))

    g = f['vector']
    V = g['block0_values'][:]
    f.close()
    return V

