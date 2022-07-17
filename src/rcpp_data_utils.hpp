#pragma once

// ------- function declaration

#include <Rcpp.h>
#include <vector>

#include "fastde_types.hpp"



//' round trip testing for dgCMatrix
//'
//' @rdname rttest_dgCMatrix
//' @param mat an R sparse compressed column matrix
//' @return R sparse compressed column matrix
//' @name rttest_dgCMatrix
//' @export
//[[Rcpp::export]]
Rcpp::dgCMatrix rttest_dgCMatrix(Rcpp::dgCMatrix const & mat);

//' round trip testing for dgCMatrix
//'
//' @rdname rttest_dgCMatrix64
//' @param mat an R sparse compressed column matrix
//' @return R sparse compressed column matrix
//' @name rttest_dgCMatrix64
//' @export
//[[Rcpp::export]]
Rcpp::dgCMatrix64 rttest_dgCMatrix64(Rcpp::dgCMatrix64 const & mat);




Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::LogicalMatrix const & _matrix, std::vector<bool> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector copy_rmatrix_to_cppvector(Rcpp::NumericMatrix const & _matrix, std::vector<double> & mat,
    size_t & nrow, size_t & ncol, size_t & nelem);

size_t copy_rvector_to_cppvector(Rcpp::LogicalVector const & _vector, std::vector<bool> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::NumericVector const & _vector, std::vector<double> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::IntegerVector const & _vector, std::vector<int> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);

size_t copy_rvector_to_cppvector(Rcpp::NumericVector const & _vector, std::vector<long> & vec, 
    size_t const & length = std::numeric_limits<size_t>::max(), size_t const & offset = 0);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix const & matrix, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);

Rcpp::StringVector copy_rsparsematrix_to_cppvectors(Rcpp::dgCMatrix64 const & matrix, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector const & _x, Rcpp::IntegerVector const & _i, Rcpp::IntegerVector const & _p, 
    std::vector<double> & x,
    std::vector<int> & i,
    std::vector<int> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


Rcpp::StringVector copy_rsparsematrix_to_cppvectors(
    Rcpp::NumericVector const & _x, Rcpp::NumericVector const & _i, Rcpp::NumericVector const & _p, 
    std::vector<double> & x,
    std::vector<long> & i,
    std::vector<long> & p,
    size_t & nrow, size_t & ncol, size_t & nelem);


// //------ these may not be needed - Rcpp can wrap?
// void import_r_common_params(SEXP as_dataframe, SEXP threads,
//     bool & _as_dataframe, int & nthreads);
// void import_de_common_params(SEXP rtype,
//     SEXP bool_param, int & type, bool & bool_val);
// void import_fc_common_params(SEXP calc_percents, 
//     SEXP min_threshold, SEXP use_expm1,
//     SEXP use_log, SEXP log_base, SEXP use_pseudocount,
//     bool & perc, double & min_thresh, bool & _use_expm1, 
//     bool & _use_log, double & _log_base, bool & _use_pseudocount);

// void import_filterfc_common_params(SEXP min_pct, 
//     SEXP min_diff_pct, SEXP logfc_threshold, SEXP only_pos, 
//     double & _min_pct, double & _min_diff_pct, double & _logfc_thresh, bool & _only_pos);
// //-----------------

Rcpp::StringVector populate_feature_names(Rcpp::StringVector const & features,
    size_t const & nfeatures);

Rcpp::NumericMatrix export_de_to_r_matrix(
    std::vector<double> const & pv,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);

Rcpp::DataFrame export_de_to_r_dataframe(
    std::vector<double> const & pv, std::string const & name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);

Rcpp::DataFrame export_fc_to_r_dataframe(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);
Rcpp::List export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<double> const & p1, std::string const & p1name,
    std::vector<double> const & p2, std::string const & p2name,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);
Rcpp::List export_fc_to_r_matrix(
    std::vector<double> const & fc, std::string const & fcname,
    std::vector<std::pair<int, size_t> > const & sorted_labels,
    Rcpp::StringVector const & features
);


class SEXP_Mat {
    protected:
        SEXP mat;
        long nrow;
        long ncol;
        int proc_count;

        SEXP copy_str_vec(SEXP src, SEXP dest, long const & count) {
            // deep copy the strings
            if (TYPEOF(src) == NILSXP) return R_NilValue;
            if (TYPEOF(dest) == NILSXP) {
                PROTECT(dest = Rf_allocVector(STRSXP, count));
                ++proc_count;
            }
            SEXP * in_name_ptr = STRING_PTR(src);
            for (long i = 0; i < count; ++i, ++in_name_ptr) {
                SET_STRING_ELT(dest, i, Rf_duplicate(*in_name_ptr));
            }
            return dest;
        }

    public:
        SEXP_Mat(long const & _nrow, long const & _ncol) :
            nrow(_nrow), ncol(_ncol), proc_count(0) {
            // https://www.r-bloggers.com/2021/04/constructing-a-sparse-matrix-class-in-rcpp/

            PROTECT(mat = Rf_allocMatrix(REALSXP, nrow, ncol));

            // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
            SEXP dimnms, rownms, colnms;
            PROTECT(dimnms = Rf_allocVector(VECSXP, 2));
            // PROTECT(rownms = Rf_allocVector(STRSXP, 0));
            // PROTECT(colnms = Rf_allocVector(STRSXP, 0));
            // SET_VECTOR_ELT(dimnms, 0, rownms);
            // SET_VECTOR_ELT(dimnms, 1, colnms);
            //  NilValue:  https://colinfay.me/r-internals/r-internal-structures.html
            SET_VECTOR_ELT(dimnms, 0, R_NilValue);
            SET_VECTOR_ELT(dimnms, 1, R_NilValue);

            Rf_setAttrib(mat, R_DimNamesSymbol, dimnms);
            proc_count = 2;
        }
        ~SEXP_Mat() { 
            if (proc_count > 0) UNPROTECT(proc_count);
            proc_count = 0;
        }

        // accessors
        // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
        inline long get_nrow() { return nrow; };
        inline long get_ncol() { return ncol; };
        inline SEXP get_rownames() { return VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 0); };
        inline SEXP get_colnames() { return VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 1); };
        inline void set_rownames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_rownames(), get_nrow());
            SET_VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 0, des);
        };
        inline void set_colnames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_colnames(), get_ncol());
            SET_VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 1, des);
        };
        inline double* get_data() { return REAL(mat); }
        inline SEXP get_sexp() { return mat; };
};


template <typename SpMAT>
class S4_SpMat;

// https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
// ======= decompose the input matrix in CSC format, S4 object with slots:
// i :  int, row numbers, 0-based.
// p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
// x :  numeric, values
// Dim:  int, 2D, sizes of full matrix
// Dimnames:  2D, names.
// factors:  ignore.

template <>
class S4_SpMat<Rcpp::dgCMatrix> {
    protected:
        Rcpp::S4 obj;
        int proc_count;
        using IDX = typename Rcpp::dgCMatrix::INDEX;

        SEXP copy_str_vec(SEXP src, SEXP dest, IDX const & count) {
            // deep copy the strings
            if (TYPEOF(src) == NILSXP) return R_NilValue;
            if (TYPEOF(dest) == NILSXP) {
                PROTECT(dest = Rf_allocVector(STRSXP, count));
                ++proc_count;
            }
            SEXP * in_name_ptr = STRING_PTR(src);
            for (IDX i = 0; i < count; ++i, ++in_name_ptr) {
                SET_STRING_ELT(dest, i, Rf_duplicate(*in_name_ptr));
            }
            return dest;
        }

    public:
        S4_SpMat(Rcpp::S4 const & _obj) : obj(_obj), proc_count(0) {};
        S4_SpMat(IDX const & nrow, IDX const & ncol, IDX const & nelem) :
            obj("dgCMatrix"), proc_count(0) {
            // https://www.r-bloggers.com/2021/04/constructing-a-sparse-matrix-class-in-rcpp/
            // create a new S4 object that matches dgCMatrix
            //  ---- create i: same size as x.
            //  ---- create x
            //  ---- p: nrow+1
            //  ---- dim:  swap original
            //  ---- dimnames: deepcopy and swap original row/col.
            SEXP i, p, x, dim, dimnms, rownms, colnms;
            PROTECT(x = Rf_allocVector(REALSXP, nelem));
            PROTECT(i = Rf_allocVector(INTSXP, nelem));
            PROTECT(p = Rf_allocVector(INTSXP, ncol + 1));
            PROTECT(dim = Rf_allocVector(INTSXP, 2));
            INTEGER(dim)[0] = nrow;
            INTEGER(dim)[1] = ncol;

            memset(INTEGER(p), 0, ncol * sizeof(int));
            INTEGER(p)[ncol] = nelem;

// https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
            PROTECT(dimnms = Rf_allocVector(VECSXP, 2));
            // PROTECT(rownms = Rf_allocVector(STRSXP, 0));
            // PROTECT(colnms = Rf_allocVector(STRSXP, 0));
            // SET_VECTOR_ELT(dimnms, 0, rownms);
            // SET_VECTOR_ELT(dimnms, 1, colnms);
            //  NilValue:  https://colinfay.me/r-internals/r-internal-structures.html
            SET_VECTOR_ELT(dimnms, 0, R_NilValue);
            SET_VECTOR_ELT(dimnms, 1, R_NilValue);
            proc_count = 5;

            obj.slot("x") = x;
            obj.slot("i") = i;
            obj.slot("p") = p;
            obj.slot("Dim") = dim;
            obj.slot("Dimnames") = dimnms;
        }
        ~S4_SpMat() { 
            if (proc_count > 0) UNPROTECT(proc_count);
            proc_count = 0;
        }

        // accessors
        // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
        inline IDX get_nelem() { return get_rowpointers()[get_ncol()]; };
        inline IDX get_nrow() { return INTEGER(obj.slot("Dim"))[0]; };
        inline IDX get_ncol() { return INTEGER(obj.slot("Dim"))[1]; };
        inline SEXP get_rownames() { return VECTOR_ELT(obj.slot("Dimnames"), 0); };
        inline SEXP get_colnames() { return VECTOR_ELT(obj.slot("Dimnames"), 1); };
        inline void set_rownames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_rownames(), get_nrow());
            SET_VECTOR_ELT(obj.slot("Dimnames"), 0, des);
        };
        inline void set_colnames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_colnames(), get_ncol());
            SET_VECTOR_ELT(obj.slot("Dimnames"), 1, des);
        };
        inline double* get_entries() { return REAL(obj.slot("x")); }
        inline IDX* get_colindices() { return INTEGER(obj.slot("i")); }
        inline IDX* get_rowpointers() { return INTEGER(obj.slot("p")); }
        inline Rcpp::S4 get_s4() { return obj; };
};


template <>
class S4_SpMat<Rcpp::dgCMatrix64> {

    protected:
        Rcpp::S4 obj;
        int proc_count;
        using IDX = typename Rcpp::dgCMatrix64::INDEX;

        SEXP copy_str_vec(SEXP src, SEXP dest, IDX const & count) {
            // deep copy the strings
            if (TYPEOF(src) == NILSXP) return R_NilValue;
            if (TYPEOF(dest) == NILSXP) {
                PROTECT(dest = Rf_allocVector(STRSXP, count));
                ++proc_count;
            }
            SEXP * in_name_ptr = STRING_PTR(src);
            for (IDX i = 0; i < count; ++i, ++in_name_ptr) {
                SET_STRING_ELT(dest, i, Rf_duplicate(*in_name_ptr));
            }
            return dest;
        }

    public:
        S4_SpMat(Rcpp::S4 const & _obj) : obj(_obj), proc_count(0) {};
        S4_SpMat(IDX const & nrow, IDX const & ncol, IDX const & nelem) :
            obj("dgCMatrix64"), proc_count(0) {
            // https://www.r-bloggers.com/2021/04/constructing-a-sparse-matrix-class-in-rcpp/
            // create a new S4 object that matches dgCMatrix
            //  ---- create i: same size as x.
            //  ---- create x
            //  ---- p: nrow+1
            //  ---- dim:  swap original
            //  ---- dimnames: deepcopy and swap original row/col.
            SEXP i, p, x, dimension, dim, dimnms, rownms, colnms;
            PROTECT(x = Rf_allocVector(REALSXP, nelem));
            PROTECT(i = Rf_allocVector(REALSXP, nelem));
            PROTECT(p = Rf_allocVector(REALSXP, ncol + 1));
            PROTECT(dimension = Rf_allocVector(REALSXP, 2));
            REAL(dimension)[0] = nrow;
            REAL(dimension)[1] = ncol;
            PROTECT(dim = Rf_allocVector(INTSXP, 2));
            INTEGER(dim)[0] = (nrow > 2147483647) ? -1 : nrow;
            INTEGER(dim)[1] = (ncol > 2147483647) ? -1 : ncol;

            memset(REAL(p), 0, ncol * sizeof(double));
            REAL(p)[ncol] = nelem;

    // https://stackoverflow.com/questions/5709940/r-extension-in-c-setting-matrix-row-column-names
            PROTECT(dimnms = Rf_allocVector(VECSXP, 2));
            // PROTECT(rownms = Rf_allocVector(STRSXP, 0));
            // PROTECT(colnms = Rf_allocVector(STRSXP, 0));
            // SET_VECTOR_ELT(dimnms, 0, rownms);
            // SET_VECTOR_ELT(dimnms, 1, colnms);
            //  NilValue:  https://colinfay.me/r-internals/r-internal-structures.html
            SET_VECTOR_ELT(dimnms, 0, R_NilValue);
            SET_VECTOR_ELT(dimnms, 1, R_NilValue);
            proc_count = 6;

            obj.slot("x") = x;
            obj.slot("i") = i;
            obj.slot("p") = p;
            obj.slot("dimension") = dimension;
            obj.slot("Dim") = dim;
            obj.slot("Dimnames") = dimnms;
        }
        ~S4_SpMat() { 
            if (proc_count > 0) UNPROTECT(proc_count);
            proc_count = 0;
        }

        // accessors
        // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
        inline IDX get_nelem() { return static_cast<IDX>(get_rowpointers()[get_ncol()]); };
        inline IDX get_nrow() { return static_cast<IDX>(REAL(obj.slot("dimension"))[0]); };
        inline IDX get_ncol() { return static_cast<IDX>(REAL(obj.slot("dimension"))[1]); };
        inline SEXP get_rownames() { return VECTOR_ELT(obj.slot("Dimnames"), 0); };
        inline SEXP get_colnames() { return VECTOR_ELT(obj.slot("Dimnames"), 1); };
        inline void set_rownames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_rownames(), get_nrow());
            SET_VECTOR_ELT(obj.slot("Dimnames"), 0, des);
        };
        inline void set_colnames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_colnames(), get_ncol());
            SET_VECTOR_ELT(obj.slot("Dimnames"), 1, des);
        };
        inline double* get_entries() { return REAL(obj.slot("x")); }
        inline double* get_colindices() { return REAL(obj.slot("i")); }
        inline double* get_rowpointers() { return REAL(obj.slot("p")); }
        inline Rcpp::S4 get_s4() { return obj; };
};
