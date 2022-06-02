#pragma once

// ==== function declaration
// from https://stackoverflow.com/questions/2808398/easily-measure-elapsed-time
#include <chrono>

template <
    class clock_t    = std::chrono::steady_clock,
    class duration_t = std::chrono::duration<double>
>
duration_t since(std::chrono::time_point<clock_t, duration_t> const& start);
