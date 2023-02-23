#pragma once 

// pure c/c++ versions

// Todo: remove all cpp11 references.  actually, using ITER makes it slower?  move this code to fastde core for now.

#include <vector>

template <typename OIT, typename XIT, typename PIT>
extern void csc_log_normalize_iter(XIT x, PIT p, size_t const & cols, double const & scale_factor, OIT out, int const & threads);

template <typename OVEC, typename XVEC, typename PVEC>
extern void csc_log_normalize_vec(XVEC const & x, PVEC const & p, size_t const & cols, double const & scale_factor, OVEC & out, int const & threads);


template <typename OIT, typename XIT, typename IIT, typename PIT>
extern void csc_clr_rows_iter(XIT x, IIT i, PIT p, size_t const & rows, size_t const & cols, OIT out, int const & threads);

template <typename OVEC, typename XVEC, typename IVEC, typename PVEC>
extern void csc_clr_rows_vec(XVEC const & x, IVEC const & i, PVEC const & p, size_t const & rows, size_t const & cols, OVEC & out, int const & threads);


template <typename OIT, typename XIT, typename IIT, typename PIT>
extern void csc_clr_cols_iter(XIT x, IIT i, PIT p, size_t const & rows, size_t const & cols, OIT out, int const & threads);

template <typename OVEC, typename XVEC, typename IVEC, typename PVEC>
extern void csc_clr_cols_vec(XVEC const & x, IVEC const & i, PVEC const & p, size_t const & rows, size_t const & cols, OVEC & out, int const & threads);



template <typename OIT, typename XIT, typename PIT>
extern void csc_relative_count_iter(XIT x, PIT p, size_t const & cols, double const & scale_factor, OIT out, int const & threads);

template <typename OVEC, typename XVEC, typename PVEC>
extern void csc_relative_count_vec(XVEC const & x, PVEC const & p, size_t const & cols, double const & scale_factor, OVEC & out, int const & threads);



