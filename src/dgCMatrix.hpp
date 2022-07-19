#pragma once

// from https://gallery.rcpp.org/articles/sparse-matrix-class/
namespace Rcpp {



// https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
// ======= decompose the input matrix in CSC format, S4 object with slots:
// i :  int, row numbers, 0-based.
// p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
// x :  numeric, values
// Dim:  int, 2D, sizes of full matrix
// Dimnames:  2D, names.
// factors:  ignore.

class dgCMatrix {
    public:
        using IDX = int;

    protected:
        S4 obj;
        int proc_count;

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
        dgCMatrix(S4 const & _obj) : obj(_obj), proc_count(0) {};
        dgCMatrix(int const & nrow, int const & ncol, IDX const & nelem) :
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
        ~dgCMatrix() { 
            if (proc_count > 0) UNPROTECT(proc_count);
            proc_count = 0;
        }

        // accessors
        // extract https://stackoverflow.com/questions/29477621/iterate-over-s4-object-slots-rcpp
        inline IDX get_nelem() const { return get_rowpointers()[get_ncol()]; };
        inline int get_nrow() const { return INTEGER(obj.slot("Dim"))[0]; };
        inline int get_ncol() const { return INTEGER(obj.slot("Dim"))[1]; };
        inline SEXP get_rownames() const { return VECTOR_ELT(obj.slot("Dimnames"), 0); };
        inline SEXP get_colnames() const { return VECTOR_ELT(obj.slot("Dimnames"), 1); };
        inline void set_rownames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_rownames(), get_nrow());
            SET_VECTOR_ELT(obj.slot("Dimnames"), 0, des);
        };
        inline void set_colnames(SEXP src) { 
            SEXP des = copy_str_vec(src, get_colnames(), get_ncol());
            SET_VECTOR_ELT(obj.slot("Dimnames"), 1, des);
        };
        inline double* get_entries() const { return REAL(obj.slot("x")); }
        inline IDX* get_colindices() const { return INTEGER(obj.slot("i")); }
        inline IDX* get_rowpointers() const { return INTEGER(obj.slot("p")); }
        inline SEXP get_entries_SEXP() const { return obj.slot("x"); }
        inline SEXP get_colindices_SEXP() const { return obj.slot("i"); }
        inline SEXP get_rowpointers_SEXP() const { return obj.slot("p"); }
        inline S4 get_S4() const { return obj; };
};

    // specialization of as
    template <> dgCMatrix as(SEXP mat);

    // specialization of wrap
    template <> SEXP wrap(const dgCMatrix& sm);

}
