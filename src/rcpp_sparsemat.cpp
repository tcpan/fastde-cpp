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



//' R Sparse Matrix Transpose
//'
//' This implementation directly constructs the new sparse matrix.  
//'     Output is also column-major, so there is random memory writes.
//' 
//' @rdname sp_transpose
//' @param sp_matrix a sparse matrix, of the form dgCMatrix
//' @return sp_matrix, transposed.
//' @name sp_transpose
//' @export
// [[Rcpp::export]]
extern SEXP sp_transpose(SEXP sp_matrix) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
    S4 obj(sp_matrix);
    SEXP i = obj.slot("i");
    SEXP p = obj.slot("p");  // ncol + 1
    SEXP x = obj.slot("x");
    SEXP dim = obj.slot("Dim");
    SEXP dimnames = obj.slot("Dimnames");
    SEXP rownms = VECTOR_ELT(dimnames, 0);
    SEXP colnms = VECTOR_ELT(dimnames, 1);
    
    
    size_t nrow = INTEGER(dim)[0];
    size_t ncol = INTEGER(dim)[1];
    size_t nelem = INTEGER(p)[ncol];   // since p is offsets, the ncol+1 entry has the total count.

    // ======= create new output
    // https://www.r-bloggers.com/2021/04/constructing-a-sparse-matrix-class-in-rcpp/
    int proc_count = 0;
    //  ---- create i: same size as x.
    //  ---- create x
    //  ---- p: nrow+1
    //  ---- dim:  swap original
    //  ---- dimnames: deepcopy and swap original row/col.
    SEXP ti, tp, tx, tdim, tdimnms, trownms, tcolnms, tmat;
    PROTECT(ti = Rf_allocVector(INTSXP, nelem));
    PROTECT(tp = Rf_allocVector(INTSXP, nrow + 1));
    PROTECT(tx = Rf_allocVector(REALSXP, nelem));
    PROTECT(tdim = Rf_allocVector(INTSXP, 2));
    PROTECT(tdimnms = Rf_allocVector(VECSXP, 2));
    PROTECT(trownms = Rf_allocVector(STRSXP, ncol));
    PROTECT(tcolnms = Rf_allocVector(STRSXP, nrow));
    proc_count += 7;
    
    // ======= do the transpose.
    // set the dimensions
    INTEGER(tdim)[0] = ncol;
    INTEGER(tdim)[1] = nrow;
    // deep copy the strings
    SEXP * in_name_ptr = STRING_PTR(rownms);
    for (size_t r = 0; r < nrow; ++r) {
        SET_STRING_ELT(tcolnms, r, Rf_duplicate(*in_name_ptr));
        ++in_name_ptr;
    }
    in_name_ptr = STRING_PTR(colnms);
    for (size_t c = 0; c < ncol; ++c) {
        SET_STRING_ELT(trownms, c, Rf_duplicate(*in_name_ptr));
        ++in_name_ptr;
    }
    SET_VECTOR_ELT(tdimnms, 0, trownms);
    SET_VECTOR_ELT(tdimnms, 1, tcolnms);
    
    // do the swap.  do random memory access instead of sorting.
    // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
    // 2. compute exclusive prefix sum as offsets in x.
    // 3. use p to get range of elements in x belonging to the same column, scatter to new x
    //     and increment offset in p.
    // 4. shift p to the right by 1 element, and set p[0] = 0
    // step 1
    memset(INTEGER(tp), 0, (nrow + 1) * sizeof(int));
    int * tp_ptr = INTEGER(tp) + 1; // this is count
    int * i_end = INTEGER(i) + nelem;  // remember i is 0-based.
    for (int * l = INTEGER(i); l != i_end; ++l) {
        ++tp_ptr[*l];
    }
    // step 2 - create max offset + 1 for each transposed row ( == new column)
    tp_ptr = INTEGER(tp);
    for (size_t off = 1; off <= nrow; ++off) {
        tp_ptr[off] += tp_ptr[off - 1];
    }
    // step 3
    int *iptr = INTEGER(i);  // 0-based.
    int *pptr = INTEGER(p) + 1;   // max (exclusive) pos for a column
    double *xptr = REAL(x);  
    tp_ptr = INTEGER(tp);       // position OFFSET in i and x.
    int *tiptr = INTEGER(ti);  // also 0-based, will contain current col ids.
    double *txptr = REAL(tx);
    int rid, cid = 0;   // column id needs to start with 0.  row ids start with 0
    int pos;
    double val;
    for (size_t l = 0; l < nelem; ++l) {
        rid = iptr[l];   // current row id (starts with 0)
        val = xptr[l];   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; l >= *pptr; ++cid, ++pptr);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = tp_ptr[rid];  // where to insert the data
        txptr[pos] = val;  // place the data
        tiptr[pos] = cid;  // place the row id (original col id. 0-based)
        ++tp_ptr[rid];  // update the offset - 1 space consumed.
    }
    // step 4
    tp_ptr = INTEGER(tp);
    memmove(tp_ptr + 1, tp_ptr, nrow * sizeof(int));
    *tp_ptr = 0;

    // ======= return
    S4 tsp_mat(std::string("dgCMatrix"));
    tsp_mat.slot("i") = ti;
    tsp_mat.slot("p") = tp;
    tsp_mat.slot("x") = tx;
    tsp_mat.slot("Dim") = tdim;
    tsp_mat.slot("Dimnames") = tdimnms;

    UNPROTECT(proc_count);
    return tsp_mat;
}
