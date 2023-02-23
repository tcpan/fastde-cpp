#pragma once

#include "normalize.hpp"


#include <algorithm>
#include <cstring>
#include <cmath>
#include <omp.h>

#include "fastde/benchmark_utils.hpp"

template <typename OIT, typename XIT, typename PIT>
extern void csc_log_normalize_iter(XIT x, PIT p, size_t const & cols, double const & scale_factor, OIT out, int const & threads) {
    
    // 2. for each element:  compute log1p(el / colsum[k] * scale_factor) = log(1.0 + el / colsum[k] * scale_factor)
    //  note that this is not the same as log1p(el) - log1p(colsum[k]) + log1p(scale_factor)
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        double colsum;

        // compute per column.
        for (auto c = offset; c < end; ++c) {
            colsum = 0;
            
            auto pstart = *(p + c);
            auto pend = *(p + c + 1);
            auto xend = x + pend;

            // compute column sum
            for (auto xptr = x + pstart; xptr != xend; ++xptr) {
                colsum += static_cast<double>(*xptr);
            }

            double factor = scale_factor / colsum;

            // compute log norm
            auto optr = out + pstart;
            for (auto xptr = x + pstart; xptr != xend; ++xptr, ++optr) {
                *optr = log1p(static_cast<double>(*xptr) * factor);
            }
        }
}


}

template <typename OVEC, typename XVEC, typename PVEC>
extern void csc_log_normalize_vec(XVEC const & x, PVEC const & p, size_t const & cols, double const & scale_factor, OVEC & out, int const & threads) {

    
    // 2. for each element:  compute log1p(el / colsum[k] * scale_factor) = log(1.0 + el / colsum[k] * scale_factor)
    //  note that this is not the same as log1p(el) - log1p(colsum[k]) + log1p(scale_factor)
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        double colsum;

        // compute per column.
        for (auto c = offset; c < end; ++c) {
            colsum = 0;
            
            auto pstart = p[c];
            auto pend = p[c+1];

            // compute column sum
            for (auto pp = pstart; pp < pend; ++pp) {
                colsum += static_cast<double>(x[pp]);
            }

            double factor = scale_factor / colsum;

            // compute log norm
            for (auto pp = pstart; pp < pend; ++pp) {
                out[pp] = log1p(static_cast<double>(x[pp]) * factor);
            }
        }
}

}



template <typename OIT, typename XIT, typename IIT, typename PIT>
extern void csc_clr_cols_iter(XIT x, IIT i, PIT p, size_t const & rows, size_t const & cols, OIT out, int const & threads) {
    
    // 2. for each element:  compute log1p(x / ( exp( sum(log1p(x[x > 0]), na.rm = TRUE) / length(x) ) ) )
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        double sumlog;   // sum of logs

        // compute per column.
        for (auto c = offset; c < end; ++c) {
            sumlog = 0;
            
            auto pstart = *(p + c);
            auto pend = *(p + c + 1);
            auto xend = x + pend;

            // compute column sum
            for (auto xptr = x + pstart; xptr != xend; ++xptr) {
                sumlog += log1p(static_cast<double>(*xptr));
            }

            double factor = 1.0 / exp(sumlog / static_cast<double>(rows));

            // compute log norm
            auto optr = out + pstart;
            for (auto xptr = x + pstart; xptr != xend; ++xptr, ++optr) {
                *optr = log1p(static_cast<double>(*xptr) * factor);
            }
        }
}

}

template <typename OVEC, typename XVEC, typename IVEC, typename PVEC>
extern void csc_clr_cols_vec(XVEC const & x, IVEC const & i, PVEC const & p, size_t const & rows, size_t const & cols, OVEC & out, int const & threads) {

    
    // 2. for each element:  compute log1p(x / ( exp( sum(log1p(x[x > 0]), na.rm = TRUE) / length(x) ) ) )
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        double sumlog;

        // compute per column.
        for (auto c = offset; c < end; ++c) {
            sumlog = 0;
            
            auto pstart = p[c];
            auto pend = p[c+1];

            // compute column sum
            for (auto pp = pstart; pp < pend; ++pp) {
                sumlog += log1p(static_cast<double>(x[pp]));
            }

            double factor = 1.0 / exp(sumlog / static_cast<double>(rows));

            // compute log norm
            for (auto pp = pstart; pp < pend; ++pp) {
                out[pp] = log1p(static_cast<double>(x[pp]) * factor);
            }
        }
}
}


template <typename OIT, typename XIT, typename IIT, typename PIT>
extern void csc_clr_rows_iter(XIT x, IIT i, PIT p, size_t const & rows, size_t const & cols, OIT out, int const & threads) {
    
    std::vector<double> rsumlog(rows, 0);

    if (threads <= 1) {
        auto nz = *(p + cols);
        auto iptr = i;
        auto xptr = x;
        // compute partial row log sums.  partitioned by columns.
        for (auto c = 0; c < nz; ++c, ++iptr, ++xptr) {
         
            rsumlog[*iptr] += log1p(static_cast<double>(*xptr));
        }

        // compute the factors
        for (auto c = 0; c < rows; ++c) {
            rsumlog[c] = 1.0 / exp( rsumlog[c] / static_cast<double>(cols));
        }

        // compute the output
        iptr = i;
        xptr = x;
        auto optr = out;

        for (auto c = 0; c < nz; ++c, ++iptr, ++xptr, ++optr) {
            *optr = log1p(*xptr * rsumlog[*iptr]);
        }
    } else {

    // 2. for each element:  compute log1p(x / ( exp( sum(log1p(x[x > 0]), na.rm = TRUE) / length(x) ) ) )
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        std::vector<double> sumlog(rows, 0);   // local sum of logs

        // compute partial row log sums.  partitioned by columns.
        for (auto c = offset; c < end; ++c) {
            
            auto pstart = *(p + c);
            auto pend = *(p + c + 1);
            auto xend = x + pend;
            auto iptr = i + pstart;

            // compute column sum
            for (auto xptr = x + pstart; xptr != xend; ++xptr, ++iptr) {
                sumlog[*iptr] += log1p(static_cast<double>(*xptr));
            }
        }
// wait for all partial sums to be completed.
#pragma omp barrier

        // now merge, partition rows by thread, then work on a different block at a time with barrier in between.
        size_t rblock = rows / threads;
        size_t rrem = rows - threads * block;
        for (auto t = 0; t < threads; ++t) {
            // {threads} number of blocks.  rotate through all blocks.  One thread per block at a time.
            int rtid = (tid + t) >= threads ? (tid + t - threads) : (tid + t);
            int rnid = rtid + 1;
            size_t roffset = rtid * rblock + (static_cast<size_t>(rtid) > rrem ? rrem : rtid);
            size_t rend = rnid * rblock + (static_cast<size_t>(rnid) > rrem ? rrem : rnid);

            // add to rsumlog
            for (auto r = roffset; r < rend; ++r) {
                rsumlog[r] += sumlog[r];
            }

// wait for all aggregations to be completed.
#pragma omp barrier
        }

        // precompute the factors
        size_t roffset = tid * rblock + (static_cast<size_t>(tid) > rrem ? rrem : tid);
        size_t rend = nid * rblock + (static_cast<size_t>(nid) > rrem ? rrem : nid);
        for (auto r = roffset; r < rend; ++r) {
            rsumlog[r] = 1.0 / exp(rsumlog[r] / static_cast<double>(cols));
        }
// wait for factors to be computed.
#pragma omp barrier

// now compute by column again.  
        for (auto c = offset; c < end; ++c) {
            // compute log norm
            auto pstart = *(p + c);
            auto pend = *(p + c + 1);
            auto xend = x + pend;
            auto optr = out + pstart;

            auto iptr = i + pstart;
            auto r = *iptr;
            for (auto xptr = x + pstart; xptr != xend; ++xptr, ++iptr, ++optr) {
                r = *iptr;
                *optr = log1p(static_cast<double>(*xptr) * rsumlog[r]);
            }
        }
}
    }

}

template <typename OVEC, typename XVEC, typename IVEC, typename PVEC>
extern void csc_clr_rows_vec(XVEC const & x, IVEC const & i, PVEC const & p, size_t const & rows, size_t const & cols, OVEC & out, int const & threads) {

    std::vector<double> rsumlog(rows, 0);

    if (threads <= 1) {
        auto nz = p[cols];

        // compute partial row log sums.  partitioned by columns.
        for (auto c = 0; c < nz; ++c) {                       
            rsumlog[i[c]] += log1p(static_cast<double>(x[c]));
        }

        // compute the factors
        for (auto c = 0; c < rows; ++c) {
            rsumlog[c] = 1.0 / exp( rsumlog[c] / static_cast<double>(cols));
        }

        // compute the output
        for (auto c = 0; c < nz; ++c) {
            out[c] = log1p(x[c] * rsumlog[i[c]]);
        }

    } else {
    // 2. for each element:  compute log1p(x / ( exp( sum(log1p(x[x > 0]), na.rm = TRUE) / length(x) ) ) )
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        std::vector<double> sumlog(rows, 0);   // local sum of logs

        // compute partial row log sums.  partitioned by columns.
        for (auto c = offset; c < end; ++c) {
                        
            auto pstart = p[c];
            auto pend = p[c+1];
            // compute column sum
            for (auto pp = pstart; pp < pend; ++pp) {
                sumlog[i[pp]] += log1p(static_cast<double>(x[pp]));
            }
        }
// wait for all partial sums to be completed.
#pragma omp barrier

        // now merge, partition rows by thread, then work on a different block at a time with barrier in between.
        size_t rblock = rows / threads;
        size_t rrem = rows - threads * block;
        for (auto t = 0; t < threads; ++t) {
            // {threads} number of blocks.  rotate through all blocks.  One thread per block at a time.
            int rtid = (tid + t) >= threads ? (tid + t - threads) : (tid + t);
            int rnid = rtid + 1;
            size_t roffset = rtid * rblock + (static_cast<size_t>(rtid) > rrem ? rrem : rtid);
            size_t rend = rnid * rblock + (static_cast<size_t>(rnid) > rrem ? rrem : rnid);

            // add to rsumlog
            for (auto r = roffset; r < rend; ++r) {
                rsumlog[r] += sumlog[r];
            }

// wait for all aggregations to be completed.
#pragma omp barrier
        }

        // precompute the factors. partition by row
        size_t roffset = tid * rblock + (static_cast<size_t>(tid) > rrem ? rrem : tid);
        size_t rend = nid * rblock + (static_cast<size_t>(nid) > rrem ? rrem : nid);
        for (auto r = roffset; r < rend; ++r) {
            rsumlog[r] = 1.0 / exp(rsumlog[r] / static_cast<double>(cols));
        }
// wait for factors to be computed.
#pragma omp barrier

// now compute by column again.  
        for (auto c = offset; c < end; ++c) {
            auto pstart = p[c];
            auto pend = p[c+1];

            // compute log norm
            for (auto pp = pstart; pp < pend; ++pp) {
                auto r = i[pp];
                out[pp] = log1p(static_cast<double>(x[pp]) * rsumlog[r]);
            }
        }
}
}

}




template <typename OIT, typename XIT, typename PIT>
extern void csc_relative_count_iter(XIT x, PIT p, size_t const & cols, double const & scale_factor, OIT out, int const & threads) {
    
    // 2. for each element:  compute (el / colsum[k] * scale_factor)
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        double colsum;

        // compute per column.
        for (auto c = offset; c < end; ++c) {
            colsum = 0;
            
            auto pstart = *(p + c);
            auto pend = *(p + c + 1);
            auto xend = x + pend;

            // compute column sum
            for (auto xptr = x + pstart; xptr != xend; ++xptr) {
                colsum += static_cast<double>(*xptr);
            }

            double factor = scale_factor / colsum;

            // compute log norm
            auto optr = out + pstart;
            for (auto xptr = x + pstart; xptr != xend; ++xptr, ++optr) {
                *optr = static_cast<double>(*xptr) * factor;
            }
        }
}

}

template <typename OVEC, typename XVEC, typename PVEC>
extern void csc_relative_count_vec(XVEC const & x, PVEC const & p, size_t const & cols, double const & scale_factor, OVEC & out, int const & threads) {

    
    // 2. for each element:  compute (el / colsum[k] * scale_factor)
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = cols / threads;
        size_t rem = cols - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        double colsum;

        // compute per column.
        for (auto c = offset; c < end; ++c) {
            colsum = 0;
            
            auto pstart = p[c];
            auto pend = p[c+1];

            // compute column sum
            for (auto pp = pstart; pp < pend; ++pp) {
                colsum += static_cast<double>(x[pp]);
            }

            double factor = scale_factor / colsum;

            // compute log norm
            for (auto pp = pstart; pp < pend; ++pp) {
                out[pp] = static_cast<double>(x[pp]) * factor;
            }
        }
}
}
