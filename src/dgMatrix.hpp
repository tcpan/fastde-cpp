#pragma once

namespace Rcpp {


class dgMatrix {
        
    protected:
        SEXP mat;
        int nrow;
        int ncol;
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
        dgMatrix(SEXP _mat) : mat(_mat) {
            SEXP dim = Rf_getAttrib( _mat, R_DimSymbol ) ;
            nrow = INTEGER(dim)[0];
            ncol = INTEGER(dim)[1];
        }
        dgMatrix(int const & _nrow, int const & _ncol) :
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
        ~dgMatrix() { 
            if (proc_count > 0) UNPROTECT(proc_count);
            proc_count = 0;
        }

        // accessors
        // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
        inline int get_nrow() const { return nrow; };
        inline int get_ncol() const { return ncol; };
        inline SEXP get_rownames() const { return VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 0); };
        inline SEXP get_colnames() const { return VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 1); };
        inline void set_rownames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_rownames(), get_nrow());
            SET_VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 0, des);
        };
        inline void set_colnames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_colnames(), get_ncol());
            SET_VECTOR_ELT(Rf_getAttrib(mat, R_DimNamesSymbol), 1, des);
        };
        inline double* get_data() const { return REAL(mat); }
        inline SEXP get_data_SEXP() const { return mat; };
        inline SEXP get_SEXP() const { return mat; };
};
        // specialization of as
    template <> dgMatrix as(SEXP mat);

    // specialization of wrap
    template <> SEXP wrap(const dgMatrix& sm);

}