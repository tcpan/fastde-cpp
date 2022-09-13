#include "fastde/benchmark_utils.tpp"

template std::chrono::duration<double> since(std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double> > const& start);
