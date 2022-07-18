#include <Rcpp.h>
#include <R.h>
#include <Rdefines.h>
using namespace Rcpp;


/*
 * wrapper for R dgCMatrix
 *
 */

#include <vector>
#include <algorithm>

// Enable C++11
// [[Rcpp::plugins(cpp11)]]

// Enable OpenMP (excludes macOS)
// [[Rcpp::plugins(openmp)]]

#include <omp.h>
#include "rcpp_data_utils.hpp"

// TODO:  [ ] add function to sparsify dense matrix.

template <typename MAT>
extern Rcpp::S4 _sp_transpose(Rcpp::S4 obj) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // input
    S4_SpMat<MAT> in(obj);
    
    // empty output 
    S4_SpMat<MAT> out(in.get_ncol(), in.get_nrow(), in.get_nelem());
    // deep copy dimnames
    out.set_colnames(in.get_rownames());
    out.set_rownames(in.get_colnames());
    

    // ======= do the transpose.
    
    // do the swap.  do random memory access instead of sorting.
    // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
    // 2. compute exclusive prefix sum as offsets in x.
    // 3. use p to get range of elements in x belonging to the same column, scatter to new x
    //     and increment offset in p.
    // 4. shift p to the right by 1 element, and set p[0] = 0
    // step 1
    auto tp_ptr = out.get_rowpointers() + 1; // this is count
    auto iptr = in.get_colindices();
    auto i_end = iptr + in.get_nelem();  // remember i is 0-based.
    for (; iptr != i_end; ++iptr) {
        ++tp_ptr[static_cast<long>(*iptr)];
    }
    // step 2 - create max offset + 1 for each transposed row ( == new column)
    tp_ptr = out.get_rowpointers();
    auto ncol = out.get_ncol();
    decltype(ncol) cid = 1;
    for (; cid <= ncol; ++cid) {
        tp_ptr[cid] += tp_ptr[cid - 1];
    }
    // step 3
    iptr = in.get_colindices();  // 0-based.
    auto pptr = in.get_rowpointers() + 1;   // max (exclusive) pos for a column
    auto xptr = in.get_entries();  
    tp_ptr = out.get_rowpointers();       // position OFFSET in i and x.
    auto tiptr = out.get_colindices();  // also 0-based, will contain current col ids.
    auto txptr = out.get_entries();
    typename std::remove_reference<decltype(iptr[0])>::type rid;   // column id needs to start with 0.  row ids start with 0
    cid = 0;
    typename std::remove_reference<decltype(xptr[0])>::type val;
    auto nelem = in.get_nelem();
    decltype(nelem) pos;
    decltype(nelem) e = 0;
    
    for (; e < nelem; ++e) {
        rid = iptr[e];   // current row id (starts with 0)
        val = xptr[e];   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; e >= *pptr; ++cid, ++pptr);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = static_cast<decltype(nelem)>(tp_ptr[static_cast<long>(rid)]);  // where to insert the data
        txptr[pos] = val;  // place the data
        tiptr[pos] = cid;  // place the row id (original col id. 0-based)
        ++tp_ptr[static_cast<long>(rid)];  // update the offset - 1 space consumed.
    }
    // step 4
    tp_ptr = out.get_rowpointers();
    memmove(tp_ptr + 1, tp_ptr, out.get_ncol() * sizeof(decltype(rid)));
    *tp_ptr = 0;

    // ======= return
    return out.get_s4();
}


//' R Sparse Matrix Transpose
//'
//' This implementation directly constructs the new sparse matrix.  
//'     Output is also column-major, so there is random memory writes.
//' 
//' @rdname rc_sp_transpose
//' @param sp_matrix a sparse matrix, of the form dgCMatrix
//' @return sp_matrix, transposed.
//' @name rc_sp_transpose
//' @export
// [[Rcpp::export]]
extern Rcpp::S4 rc_sp_transpose(Rcpp::S4 sp_matrix) {
    return _sp_transpose<Rcpp::dgCMatrix>(sp_matrix);
}


//' R Sparse Matrix Transpose
//'
//' This implementation directly constructs the new sparse matrix.  
//'     Output is also column-major, so there is random memory writes.
//' 
//' @rdname rc_sp64_transpose
//' @param obj a sparse matrix, of the form dgCMatrix
//' @return sp_matrix, transposed.
//' @name rc_sp64_transpose
//' @export
// [[Rcpp::export]]
extern Rcpp::S4 rc_sp64_transpose(Rcpp::S4 obj) {
    return _sp_transpose<Rcpp::dgCMatrix64>(obj);

}


template <typename MAT>
extern SEXP _sp_to_dense(Rcpp::S4 obj) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
    S4_SpMat<MAT> in(obj);

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    // ======= create new output and initialize
    SEXP_Mat out(in.get_nrow(), in.get_ncol());
    out.set_rownames(in.get_rownames());
    out.set_colnames(in.get_colnames());

    
    double * data = out.get_data();
    memset(data, 0, out.get_nrow() * out.get_ncol() * sizeof(double));

    auto cp = in.get_rowpointers();
    auto ri = in.get_colindices();
    double * xptr = in.get_entries();
    // size_t r, offset;
    // size_t mx = cp[ncol];
    // size_t c = 0;
    // size_t c_end = cp[1];
    // for (size_t e = 0; e < mx; ++e ) {
    //     while (e >= c_end ) {
    //         ++c;
    //         c_end = cp[c + 1];
    //         data = REAL(dense) + c * nrow;
    //     }
    //     data[ri[e]] = xptr[e];
    // }
    long o_end = static_cast<long>(cp[0]);
    long o;
    auto ncol = in.get_ncol();
    auto nrow = in.get_nrow();
    decltype(ncol) c = 0;
    for (; c < ncol; ++c) {
        o = o_end;
        o_end = static_cast<long>(cp[c + 1]);
        data = out.get_data() + c * nrow;
        for (; o < o_end; ++o) {
            data[static_cast<long>(ri[o])] = xptr[o];
        }
    }

    return out.get_sexp();
}


//' R Sparse To Dense Matrix
//'
//' This implementation directly constructs the new dense matrix.  
//'     There is random memory writes.
//' 
//' @rdname rc_sp_to_dense
//' @param obj a sparse matrix, of the form dgCMatrix
//' @return matrix dense matrix.
//' @name rc_sp_to_dense
//' @export
// [[Rcpp::export]]
extern SEXP rc_sp_to_dense(S4 obj) {
    return _sp_to_dense<Rcpp::dgCMatrix>(obj);
}



//' R Sparse To Dense Matrix
//'
//' This implementation directly constructs the new dense matrix.  
//'     There is random memory writes.
//' 
//' @rdname rc_sp64_to_dense
//' @param obj a sparse matrix, of the form dgCMatrix
//' @return matrix dense matrix.
//' @name rc_sp64_to_dense
//' @export
// [[Rcpp::export]]
extern SEXP rc_sp64_to_dense(S4 obj) {
    return _sp_to_dense<Rcpp::dgCMatrix64>(obj);
}


template <typename MAT>
extern SEXP _sp_to_dense_transposed(Rcpp::S4 obj) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
    S4_SpMat<MAT> in(obj);

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    // ======= create new output and initialize
    SEXP_Mat out(in.get_ncol(), in.get_nrow());
    out.set_rownames(in.get_colnames());
    out.set_colnames(in.get_rownames());

    double * data = out.get_data();
    memset(data, 0, out.get_nrow() * out.get_ncol() * sizeof(double));

    auto cp = in.get_rowpointers();
    auto ri = in.get_colindices();
    double * xptr = in.get_entries();
    // size_t r, offset;
    // size_t mx = cp[ncol];
    // size_t c = 0;
    // size_t c_end = cp[1];
    // for (size_t e = 0; e < mx; ++e ) {
    //     while (e >= c_end ) {
    //         ++c;
    //         c_end = cp[c + 1];
    //         data = REAL(dense) + c * nrow;
    //     }
    //     data[ri[e]] = xptr[e];
    // }
    long o_end = static_cast<long>(cp[0]);
    long o;
    auto ncol = in.get_ncol();
    auto nrow = in.get_nrow();
    decltype(ncol) c = 0;
    for (; c < ncol; ++c) {
        o = o_end;
        o_end = static_cast<long>(cp[c + 1]);
        data = out.get_data() + c;
        for (; o < o_end; ++o) {
            data[static_cast<long>(ri[o]) * ncol] = xptr[o];
        }
    }

    return out.get_sexp();
}


//' R Sparse To Dense Matrix
//'
//' This implementation directly constructs the new dense matrix.  
//'     There is random memory writes.
//' 
//' @rdname rc_sp_to_dense_transposed
//' @param obj a sparse matrix, of the form dgCMatrix
//' @return matrix dense matrix.
//' @name rc_sp_to_dense_transposed
//' @export
// [[Rcpp::export]]
extern SEXP rc_sp_to_dense_transposed(S4 obj) {
    return _sp_to_dense_transposed<Rcpp::dgCMatrix>(obj);
}

//' R Sparse To Dense Matrix
//'
//' This implementation directly constructs the new dense matrix.  
//'     There is random memory writes.
//' 
//' @rdname rc_sp64_to_dense_transposed
//' @param obj a sparse matrix, of the form dgCMatrix
//' @return matrix dense matrix.
//' @name rc_sp64_to_dense_transposed
//' @export
// [[Rcpp::export]]
extern SEXP rc_sp64_to_dense_transposed(S4 obj) {
    return _sp_to_dense_transposed<Rcpp::dgCMatrix64>(obj);
}
