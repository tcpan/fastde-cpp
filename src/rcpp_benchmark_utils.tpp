#pragma once

// ==== function definition
#include "rcpp_benchmark_utils.hpp"


template < class clock_t, class duration_t >
duration_t since(std::chrono::time_point<clock_t, duration_t> const& start)
{
    return std::chrono::duration_cast<duration_t>(clock_t::now() - start);
}

SEXP start_profiler(SEXP str) {
  ProfilerStart(Rcpp::as<const char*>(str));
  return R_NilValue;
}

SEXP stop_profiler() {
  ProfilerStop();
  return R_NilValue;
}