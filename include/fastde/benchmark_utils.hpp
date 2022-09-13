#pragma once

// ==== function declaration
// from https://stackoverflow.com/questions/2808398/easily-measure-elapsed-time
#include <chrono>
// #include "gperftools/profiler.h"

template <
    class clock_t    = std::chrono::steady_clock,
    class duration_t = std::chrono::duration<double>
>
duration_t since(std::chrono::time_point<clock_t, duration_t> const& start);

// //' start profiler
// //'
// //' @rdname start_profiler
// //' @param str name of output file
// //' @return Nil value
// //' @name start_profiler
// //' @export
// [[cpp11::register]]
// extern SEXP start_profiler(SEXP str);

// //' stop profiler
// //'
// //' @rdname stop_profiler
// //' @name stop_profiler
// //' @export
// [[cpp11::register]]
// extern SEXP stop_profiler();