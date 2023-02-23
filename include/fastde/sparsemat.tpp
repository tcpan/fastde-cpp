#include "fastde/sparsemat.hpp"

/*
 * wrapper for R dgCMatrix
 *
 */

#include <algorithm>
#include <cstring>

#include <omp.h>

#include "fastde/benchmark_utils.hpp"


// #include <atomic>
// #include <parallel/algorithm>
// VERY SLOW
// // NOTe:  there is no formal definition of sparse matrix.
// // input is column major, so i has the row ids, and p is per column.
// template <typename XT, typename IT, typename PT, typename IT2>
// extern cpp11::writable::list csc_transpose_sort(
//     cpp11::r_vector<XT> const & x, 
//     cpp11::r_vector<IT> const & i, 
//     cpp11::r_vector<PT> const & p, IT2 const & nrow, IT2 const & ncol, int const & threads) {

//     using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;
//     // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
//     // ======= decompose the input matrix in CSC format, S4 object with slots:
//     // i :  int, row numbers, 0-based.
//     // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
//     // x :  numeric, values
//     // Dim:  int, 2D, sizes of full matrix
//     // Dimnames:  2D, names.
//     // factors:  ignore.
    
//     // input
//     size_t nelem = x.size();
//   std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;


//   start = std::chrono::steady_clock::now();

//     std::vector<std::tuple<IT, PT2, XT>> temp;
//     temp.resize(nelem);

//     // copy to temp
// #pragma omp parallel num_threads(threads)
// {   
//     int tid = omp_get_thread_num();
//     size_t block = ncol / threads;
//     size_t rem = ncol - threads * block;
//     size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
//     int nid = tid + 1;
//     size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

//     IT2 i_start = p[offset], i_end;
//     for (; offset < end; ++offset) {
//         i_end = p[offset + 1];
//         for (; i_start < i_end; ++i_start) {
//             temp[i_start] = {i[i_start], offset, x[i_start]};
//         }
//     }
// }


//     // parallel sort
//     __gnu_parallel::sort(temp.begin(), temp.end(), 
//         [](std::tuple<IT, PT2, XT> const & lhs, std::tuple<IT, PT2, XT> const & rhs){
//             return (std::get<0>(lhs) == std::get<0>(rhs)) ? (std::get<1>(lhs) < std::get<1>(rhs)) : (std::get<0>(lhs) < std::get<0>(rhs));
//     });


//   Rprintf("[TIME] sp_transpose_par 1 Elapsed(ms)= %f\n", since(start).count());

//   start = std::chrono::steady_clock::now();
//     // looking for starting offset for the rows. 
//     // need to do below in case 1 row spans multiple threads.
//     // collect from each thread, the first row encountered's offset position.
//     std::vector<std::pair<IT, PT2>> boundaries;
//     boundaries.resize(threads);
// #pragma omp parallel num_threads(threads)
// {   
//     int tid = omp_get_thread_num();
//     size_t block = nelem / threads;
//     size_t rem = nelem - threads * block;
//     size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);

//     boundaries[tid] = {std::get<0>(temp[offset]), offset};
// }


//   Rprintf("[TIME] sp_transpose_par 2 Elapsed(ms)= %f\n", since(start).count());

//   start = std::chrono::steady_clock::now();
//     // empty output 
//     cpp11::writable::r_vector<XT> tx(nelem); 
//     cpp11::writable::r_vector<IT> ti(nelem);   // as many as there are values 
//     cpp11::writable::r_vector<PT> tp(nrow + 1);       // number of rows + 1.
    
//     std::fill(tp.begin(), tp.end(), nelem);  // initialize tp so we can do min.

//     // do a linear scan of the thread boundaries for min
//     IT rid;
//     for (int t = 0; t < threads; ++t) {
//         rid = boundaries[t].first;
//         tp[rid] = (boundaries[t].second < tp[rid]) ? boundaries[t].second : tp[rid];
//     }
//     // now that we've set a minimum offsets for the boundary cases, we can blanketly compare and update -
//     // only 1 thread will have temp elements with offsets that are smaller to update tp
//     // so no thread contention.

//   Rprintf("[TIME] sp_transpose_par 3 Elapsed(ms)= %f\n", since(start).count());

//   start = std::chrono::steady_clock::now();
//     // now copy out.
// #pragma omp parallel num_threads(threads)
// {   
//     int tid = omp_get_thread_num();
//     size_t block = nelem / threads;
//     size_t rem = nelem - threads * block;
//     size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
//     int nid = tid + 1;
//     size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

//     IT rid2;
//     for (; offset < end; ++offset) {
//         tx[offset] = std::get<2>(temp[offset]);
//         ti[offset] = std::get<1>(temp[offset]);
//         rid2 = std::get<0>(temp[offset]);
//         if (offset < tp[rid2]) tp[rid2] = offset;
//     }
// }
//   Rprintf("[TIME] sp_transpose_par 4 Elapsed(ms)= %f\n", since(start).count());

//   start = std::chrono::steady_clock::now();

//     // scan tp to replace remaining "nelem" value, i.e. empty columns
//     PT v = 0;
//     for (size_t e = 0; e < nrow; ++e) {
//         if (tp[e] == nelem) tp[e] = v;
//         else v = tp[e];
//     }
//   Rprintf("[TIME] sp_transpose_par 5 Elapsed(ms)= %f\n", since(start).count());

//   start = std::chrono::steady_clock::now();
//     // ======= return
//     cpp11::named_arg _tx("x"); _tx = tx;
//     cpp11::named_arg _ti("i"); _ti = ti;
//     cpp11::named_arg _tp("p"); _tp = tp;
//     cpp11::writable::list out( { _tx, _ti, _tp} );
//   Rprintf("[TIME] sp_transpose_par 6 Elapsed(ms)= %f\n", since(start).count());

//     return out;
// }


// NOTE:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XITER, typename IITER, typename PITER, typename IT2,
    typename XITER2, typename IITER2, typename PITER2>
extern void csc_transpose_csc(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol,
    XITER2 tx,  // nelem
    IITER2 ti,  // nelem
    PITER2 tp,  // nrow+1
    int const & threads) {

    using XT = typename std::iterator_traits<XITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using iPT = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;
    // std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;

    iPT nelem = *(p + ncol);

    if (threads > 1) {

        // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
        // ======= decompose the input matrix in CSC format, S4 object with slots:
        // i :  int, row numbers, 0-based.
        // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
        // x :  numeric, values
        // Dim:  int, 2D, sizes of full matrix
        // Dimnames:  2D, names.
        // factors:  ignore.
    
        // either:  per thread summary,   this would still use less memory than sortiing the whole thing.
        // bin by row in random (thread) order, then sort per row -  this would be n log n - n log t
    
        // start = std::chrono::steady_clock::now();

        // set up per thread offsets.
        // should have
        //      r0      r1      r2      r3  ...
        // t0   0       s1      s3
        // t1   c1      s1+d1
        // t2   c1+c2   s1+d1+d2
        // ...
        // tn   s1      s2
        PT * gp = reinterpret_cast<PT *>(malloc((threads + 1) * (nrow + 1) * sizeof(PT)));
        std::vector<PT *> lps(threads + 1);
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        lps[tid] = gp + tid * (nrow + 1);
        memset(lps[tid], 0, (nrow + 1) * sizeof(PT));
}
        lps[threads] = gp + threads * (nrow + 1);
        memset(lps[threads], 0, (nrow + 1) * sizeof(PT));

        // Rprintf("[TIME] sp_transpose_par 2 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // ======= do the transpose.
        
        // do the swap.  do random memory access instead of sorting.
        // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
        // 2. compute exclusive prefix sum as offsets in x.
        // 3. use p to get range of elements in x belonging to the same column, scatter to new x
        //     and increment offset in p.
        // 4. shift p to the right by 1 element, and set p[0] = 0

        // step 1: do local count.  store in thread + 1, row r (not r+1).  partition elements.
        // i.e. compute
        //      r0      r1      r2      r3  ...
        // t0   0       0       0
        // t1   c1      d1
        // t2   c2      d2
        // ...
        // tn   cn      dn
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nelem / threads;
        size_t rem = nelem - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        auto istart = i + offset;
        auto iend = i + end;
        for (; istart != iend; ++istart) {
            ++lps[nid][static_cast<IT2>(*istart)];
        }
}
        // Rprintf("[TIME] sp_transpose_par 3 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 1.1:  for each row, prefix sum for threads, store in thread t+1, row r.  partition rows.
        // i.e. compute
        //      r0      r1      r2      r3  ...
        // t0   0       0       0
        // t1   c1      d1
        // t2   c1+c2   d1+d2
        // ...
        // tn   c(1..n) d(1..n)
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = (nrow+1) / threads;   // each thread handles a block of rows.
        size_t rem = (nrow+1) - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        for (int t = 0; t < threads; ++t) {  // linear scan, for hardware prefetching.
            for (size_t r = offset; r < end; ++r) {
                lps[t+1][r] += lps[t][r];
            }
        }
        // at the end, lps[thread] has total counts per row.
}
        // Rprintf("[TIME] sp_transpose_par 4 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 2: global prefix sum of lps[thread] by row r,  store output in lps[0], row r+1
        // linear..
        // also step 4. copy to output p array.
        // step 1.1:  for each row, prefix sum for threads, store in thread t+1, row r.  partition rows.
        // i.e. compute
        //      r0      r1      r2      r3  ...
        // t0   0       s1      s1+s2
        // t1   c1      d1
        // t2   c1+c2   d1+d2
        // ...
        // tn   s1      s2      s3
        for (IT2 r = 0; r < nrow; ++r) {
            lps[0][r+1] = lps[0][r] + lps[threads][r];
            *(tp+r) = lps[0][r];
        }
        *(tp+nrow) = lps[0][nrow];

        // Rprintf("[TIME] sp_transpose_par 5 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 2.1: add global prefix to local.  do in parallel.  each thread can do independently.
        //      r0      r1      r2      r3  ...
        // t0   0       s1      s3
        // t1   c1      s1+d1
        // t2   c1+c2   s1+d1+d2
        // ...
        // tn   s1      s1+s2
#pragma omp parallel num_threads(threads)
{   
        int nid = omp_get_thread_num() + 1;
        for (IT2 r = 0; r < nrow; ++r) {
            lps[nid][r] += lps[0][r];
        }
}
        // per thread we now have the starting offset for writing.

        // Rprintf("[TIME] sp_transpose_par 6 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 3.  use the per thread offsets to write out.
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nelem / threads;
        size_t rem = nelem - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        IT2 rid;   // column id needs to start with 0.  row ids start with 0
        XT val;
        // need to search for cid based on offset.
        auto pptr = std::upper_bound(p, p + ncol + 1, offset);
        IT2 cid = std::distance(p, pptr) - 1;
        iPT pos;
        
        auto iiter = i + offset;
        auto xiter = x + offset;
        for (; offset < end; ++offset, ++iiter, ++xiter) {
            rid = *iiter;   // current row id (starts with 0)
            val = *xiter;   // current value
            // if the current element pos reaches first elem of next column (*pptr),
            // then go to next column (increment cid and pptr).
            for (; offset >= static_cast<size_t>(*pptr); ++cid, ++pptr);  // current column id

            // now copy and update.
            // curr pos is the offset for the transposed row (new col), in tp.
            // note we are using tp array to track current offset.
            pos = lps[tid][rid];  // where to insert the data
            *(tx + pos) = val;  // place the data
            *(ti + pos) = cid;  // place the row id (original col id. 0-based)
            ++lps[tid][rid];  // update the offset - 1 space consumed.
        }
}
        free(gp);
        // Rprintf("[TIME] sp_transpose_par 7 Elapsed(ms)= %f\n", since(start).count());

    } else {  // single thread.

        // initialize
        auto tpptr = tp;
        for (size_t r = 0; r <= nrow; ++r, ++tpptr) { *tpptr = 0; }
        // std::fill(tp, tp + nrow + 1, 0);

        // do the swap.  do random memory access instead of sorting.
        // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
        // 2. compute exclusive prefix sum as offsets in x.
        // 3. use p to get range of elements in x belonging to the same column, scatter to new x
        //     and increment offset in p.
        // 4. shift p to the right by 1 element, and set p[0] = 0
        // step 1
        auto iptr = i;
        auto i_end = i + nelem;  // remember i is 0-based.
        for (; iptr != i_end; ++iptr) {
            *(tp + static_cast<IT2>(*iptr) + 1) += 1;
        }
        // step 2 - create max offset + 1 for each transposed row ( == new column)
        auto tptr = tp;
        auto t_end = tp + nrow;
        for (; tptr != t_end; ++tptr) {
            *(tptr + 1) += *tptr;
        }
        // step 3  scatter
        IT2 rid;   // column id needs to start with 0.  row ids start with 0
        auto pptr = p + 1;  // compare to end of col ptr
        IT2 cid = 0;
        XT val;
        iPT pos;
        tpptr = tp;
        for (iPT e = 0; e < nelem; ++e, ++i, ++x) {
            rid = *i;   // current row id (starts with 0)
            val = *x;   // current value
            // if the current element pos reaches first elem of next column (*pptr),
            // then go to next column (increment cid and pptr).
            for (; e >= static_cast<iPT>(*pptr); ++cid, ++pptr);  // current column id

            // now copy and update.
            // curr pos is the offset for the transposed row (new col), in tp.
            // note we are using tp array to track current offset.
            pos = *(tpptr + rid);  // where to insert the data
            *(tx + pos) = val;  // place the data
            *(ti + pos) = cid;  // place the row id (original col id. 0-based)
            *(tpptr + rid) += 1;  // update the offset - 1 space consumed.
        }
        // step 4.  shift the tp values forward by 1 position.
        // std::copy_backward(tp, tp + nrow, tp + nrow + 1);
        size_t temp = 0, temp2; 
        for (size_t r = 0; r <= nrow; ++r, ++tp) {
            temp2 = *tp;
            *tp = temp;
            temp = temp2;
        }

    }

}




// NOTE:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XVEC, typename IVEC, typename PVEC, typename IT2,
    typename XVEC2, typename IVEC2, typename PVEC2>
extern void csc_transpose_csc_vec(
    XVEC const & x, 
    IVEC const & i, 
    PVEC const & p, 
    IT2 const & nrow, IT2 const & ncol,
    XVEC2 & tx,  // nelem
    IVEC2 & ti,  // nelem
    PVEC2 & tp,  // nrow+1
    int const & threads) {


    size_t nelem = p[ncol];

    if (threads > 1) {

        // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
        // ======= decompose the input matrix in CSC format, S4 object with slots:
        // i :  int, row numbers, 0-based.
        // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
        // x :  numeric, values
        // Dim:  int, 2D, sizes of full matrix
        // Dimnames:  2D, names.
        // factors:  ignore.
    
        // either:  per thread summary,   this would still use less memory than sortiing the whole thing.
        // bin by row in random (thread) order, then sort per row -  this would be n log n - n log t
    
        // start = std::chrono::steady_clock::now();

        // set up per thread offsets.
        // should have
        //      r0      r1      r2      r3  ...
        // t0   0       s1      s3
        // t1   c1      s1+d1
        // t2   c1+c2   s1+d1+d2
        // ...
        // tn   s1      s2

        using PT = typename std::remove_cv<typename std::remove_reference<decltype(tp[0])>::type>::type;
        std::vector<std::vector<PT>> lps(threads + 1);

#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        lps[tid].resize(nrow + 1, 0);
}
        lps[threads].resize(nrow + 1, 0);

        // Rprintf("[TIME] sp_transpose_par 2 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // ======= do the transpose.
        
        // do the swap.  do random memory access instead of sorting.
        // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
        // 2. compute exclusive prefix sum as offsets in x.
        // 3. use p to get range of elements in x belonging to the same column, scatter to new x
        //     and increment offset in p.
        // 4. shift p to the right by 1 element, and set p[0] = 0

        // step 1: do local count.  store in thread + 1, row r (not r+1).  partition elements.
        // i.e. compute
        //      r0      r1      r2      r3  ...
        // t0   0       0       0
        // t1   c1      d1
        // t2   c2      d2
        // ...
        // tn   cn      dn
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nelem / threads;
        size_t rem = nelem - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        for (size_t e = offset; e < end; ++e) {
            ++lps[nid][static_cast<size_t>(i[e])];
        }
}
        // Rprintf("[TIME] sp_transpose_par 3 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 1.1:  for each row, prefix sum for threads, store in thread t+1, row r.  partition rows.
        // i.e. compute
        //      r0      r1      r2      r3  ...
        // t0   0       0       0
        // t1   c1      d1
        // t2   c1+c2   d1+d2
        // ...
        // tn   c(1..n) d(1..n)
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = (nrow+1) / threads;   // each thread handles a block of rows.
        size_t rem = (nrow+1) - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        for (int t = 0; t < threads; ++t) {  // linear scan, for hardware prefetching.
            for (size_t r = offset; r < end; ++r) {
                lps[t+1][r] += lps[t][r];
            }
        }
        // at the end, lps[thread] has total counts per row.
}
        // Rprintf("[TIME] sp_transpose_par 4 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 2: global prefix sum of lps[thread] by row r,  store output in lps[0], row r+1
        // linear..
        // also step 4. copy to output p array.
        // step 1.1:  for each row, prefix sum for threads, store in thread t+1, row r.  partition rows.
        // i.e. compute
        //      r0      r1      r2      r3  ...
        // t0   0       s1      s1+s2
        // t1   c1      d1
        // t2   c1+c2   d1+d2
        // ...
        // tn   s1      s2      s3
        for (size_t r = 0; r < nrow; ++r) {
            lps[0][r+1] = lps[0][r] + lps[threads][r];
            tp[r] = lps[0][r];
        }
        tp[nrow] = lps[0][nrow];

        // Rprintf("[TIME] sp_transpose_par 5 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 2.1: add global prefix to local.  do in parallel.  each thread can do independently.
        //      r0      r1      r2      r3  ...
        // t0   0       s1      s3
        // t1   c1      s1+d1
        // t2   c1+c2   s1+d1+d2
        // ...
        // tn   s1      s1+s2
#pragma omp parallel num_threads(threads)
{   
        int nid = omp_get_thread_num() + 1;
        for (size_t r = 0; r < nrow; ++r) {
            lps[nid][r] += lps[0][r];
        }
}
        // per thread we now have the starting offset for writing.

        // Rprintf("[TIME] sp_transpose_par 6 Elapsed(ms)= %f\n", since(start).count());

        // start = std::chrono::steady_clock::now();

        // step 3.  use the per thread offsets to write out.
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nelem / threads;
        size_t rem = nelem - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        size_t rid;   // column id needs to start with 0.  row ids start with 0
        auto val = x[0];
        // need to search for cid based on offset.
        auto pptr = std::upper_bound(p.begin(), p.end(), offset);
        size_t ppos = std::distance(p.begin(), pptr), cid;
        size_t pos;
        
        for (size_t e = offset ; e < end; ++e) {
            rid = i[e];   // current row id (starts with 0)
            val = x[e];   // current value
            // if the current element pos reaches first elem of next column (*pptr),
            // then go to next column (increment cid and pptr).
            for (; e >= static_cast<size_t>(p[ppos]); ++ppos);  // current column id
            cid = ppos - 1;

            // now copy and update.
            // curr pos is the offset for the transposed row (new col), in tp.
            // note we are using tp array to track current offset.
            pos = lps[tid][rid];  // where to insert the data
            tx[pos] = val;  // place the data
            ti[pos] = cid;  // place the row id (original col id. 0-based)
            ++lps[tid][rid];  // update the offset - 1 space consumed.
        }
}
        // Rprintf("[TIME] sp_transpose_par 7 Elapsed(ms)= %f\n", since(start).count());

    } else {  // single thread.

        // initialize
        for (size_t r = 0; r <= nrow; ++r) { tp[r] = 0; }
        // tp.resize(nrow + 1, 0);
        // std::fill(tp.begin(), tp.end(), 0);

        // do the swap.  do random memory access instead of sorting.
        // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
        // 2. compute exclusive prefix sum as offsets in x.
        // 3. use p to get range of elements in x belonging to the same column, scatter to new x
        //     and increment offset in p.
        // 4. shift p to the right by 1 element, and set p[0] = 0
        // step 1
        for (size_t e = 0; e < nelem; ++e) {
            ++tp[i[e] + 1];
        }
        // step 2 - create max offset + 1 for each transposed row ( == new column)
        for (size_t c = 0; c < nrow; ++c) {
            tp[c+1] += tp[c];
        }
        // step 3  scatter
        size_t rid;   // column id needs to start with 0.  row ids start with 0
        size_t cid = 0, ppos = 1;
        auto val = x[0];
        size_t pos;
        for (size_t e = 0; e < nelem; ++e) {
            rid = i[e];   // current row id (starts with 0)
            val = x[e];   // current value
            // if the current element pos reaches first elem of next column (*pptr),
            // then go to next column (increment cid and pptr).
            for (; e >= static_cast<size_t>(p[ppos]); ++ppos);  // current column id
            cid = ppos - 1;

            // now copy and update.
            // curr pos is the offset for the transposed row (new col), in tp.
            // note we are using tp array to track current offset.
            pos = tp[rid];  // where to insert the data
            tx[pos] = val;  // place the data
            ti[pos] = cid;  // place the row id (original col id. 0-based)
            ++tp[rid];  // update the offset - 1 space consumed.
        }
        // step 4.  shift the tp values forward by 1 position.
        // std::copy_backward(tp, tp + nrow, tp + nrow + 1);
        size_t temp = 0, temp2; 
        for (size_t r = 0; r <= nrow; ++r) {
            temp2 = tp[r];
            tp[r] = temp;
            temp = temp2;
        }

    }

}


// compressed sparse column to column major.
template <typename OITER, typename XITER, typename IITER,
    typename PITER, typename IT2>
extern void csc_to_dense_c(
    XITER x, 
    IITER i,
    PITER p, 
    IT2 const & nrow, IT2 const & ncol, 
    OITER out,
    int const & threads) {


    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    using PT = typename std::iterator_traits<PITER>::value_type;
    using iPT = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

    // ======= create new output and initialize
    // std::fill(out, out + nrow * ncol, 0);
    size_t nz = nrow * ncol;
    auto optr = out;
    for (size_t e = 0; e < nz; ++e, ++optr) {*optr = 0; }

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    if (threads == 1) {
        IT2 r;
        iPT istart, iend;
        auto pptr = p;
        auto xptr = x;
        auto iptr = i;
        auto ipend = i;
        // auto optr = out;
        for (size_t c = 0; c < ncol; ++c) {
            istart = *pptr;
            ++pptr;
            iend = *pptr;

            xptr = x + istart;
            iptr = i + istart;
            ipend = i + iend;
            OITER optr(out + c * nrow);
            for (; iptr != ipend; ++xptr, ++iptr) {
                r = *iptr;
                *(optr + r) = *xptr;
            }
        }

    } else {
    // now iterate and fill   
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        IT2 r;
        iPT istart, iend;
        auto pptr = p + offset;
        auto iptr = i;
        auto ipend = i;
        auto xptr = x;
        // auto optr = out;
        for (; offset < end; ++offset) {
            istart = *pptr;
            ++pptr;
            iend = *pptr;

            xptr = x + istart;
            iptr = i + istart;
            ipend = i + iend;
            OITER optr(out + offset * nrow);
            for (; iptr != ipend; ++xptr, ++iptr) {
                r = *iptr;
                *(optr + r) = *xptr;
            }
        }
}
    }
}


// compressed sparse column to column major.
template <typename OVEC, typename XVEC, typename IVEC,
    typename PVEC, typename IT2>
extern void csc_to_dense_c_vec(
    XVEC const & x, 
    IVEC const & i,
    PVEC const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    OVEC & out,
    int const & threads) {


    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.

    // ======= create new output and initialize
    size_t nz = nrow * ncol;
    for (size_t e = 0; e < nz; ++e) { out[e] = 0; }
    // out.clear();
    // out.resize(nrow * ncol, 0);

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 
    
    if (threads == 1) {
        size_t istart, iend = p[0];
        auto r = i[0];
        size_t off = 0;

        // auto optr = out;
        for (size_t c = 0; c < ncol; ++c) {
            istart = iend;
            iend = p[c+1];
            off = c * nrow;

            for (size_t e = istart; e < iend; ++e) {
                r = i[e];
                out[r + off] = x[e];
            }
        }

    } else {
    // now iterate and fill   
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        auto r = i[0];
        size_t istart, iend = p[offset];
        size_t off;
        // auto optr = out;
        for (size_t c = offset; c < end; ++c) {
            istart = iend;
            iend = p[c + 1];
            off = c * nrow;

            for (size_t e = istart; e < iend; ++e) {
                r = i[e];
                out[r + off] = x[e];
            }
        }
}
    }
}

template <typename OITER, typename XITER, typename IITER,
    typename PITER, typename IT2>
extern void csc_to_dense_transposed_c(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol, 
    OITER out,
    int const & threads) {


    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    using PT = typename std::iterator_traits<PITER>::value_type;
    using iPT = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

    // ======= create new output and initialize
    // std::fill(out, out + nrow * ncol, 0);
    size_t nz = nrow * ncol;
    auto optr = out;
    for (size_t e = 0; e < nz; ++e, ++optr) {*optr = 0; }

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    if (threads == 1) {
        IT2 r;
        iPT istart, iend;
        auto pptr = p;
        auto xptr = x;
        auto iptr = i;
        auto ipend = i;
        // auto optr = out;
        for (size_t c = 0; c < ncol; ++c) {
            istart = *pptr;
            ++pptr;
            iend = *pptr;

            xptr = x + istart;
            iptr = i + istart;
            ipend = i + iend;
            OITER optr(out + c);
            for (; iptr != ipend; ++xptr, ++iptr) {
                r = *iptr;
                *(optr + r * ncol) = *xptr;
            }
        }

    } else {
    // now iterate and fill   
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        IT2 r;
        iPT istart, iend;
        auto pptr = p + offset;
        auto iptr = i;
        auto ipend = i;
        auto xptr = x;
        // auto optr = out;
        for (; offset < end; ++offset) {
            istart = *pptr;
            ++pptr;
            iend = *pptr;

            xptr = x + istart;
            iptr = i + istart;
            ipend = i + iend;
            OITER optr(out + offset);
            for (; iptr != ipend; ++xptr, ++iptr) {
                r = *iptr;
                *(optr + r * ncol) = *xptr;
            }
        }
}
    }
}


template <typename OVEC, typename XVEC, typename IVEC,
    typename PVEC, typename IT2>
extern void csc_to_dense_transposed_c_vec(
    XVEC const & x, 
    IVEC const & i, 
    PVEC const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    OVEC & out,
    int const & threads) {


    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    

    // ======= create new output and initialize
    size_t nz = nrow * ncol;
    for (size_t e = 0; e < nz; ++e) { out[e] = 0; }
    // out.clear();
    // out.resize(nrow * ncol, 0);

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 
    
    if (threads == 1) {
        size_t istart, iend = p[0];
        auto r = i[0];

        // auto optr = out;
        for (size_t c = 0; c < ncol; ++c) {
            istart = iend;
            iend = p[c+1];

            for (size_t e = istart; e < iend; ++e) {
                r = i[e];
                out[c + r * ncol] = x[e];
            }
        }

    } else {
    // now iterate and fill   
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        auto r = i[0];
        size_t istart, iend = p[offset];
        // auto optr = out;
        for (size_t c = offset; c < end; ++c) {
            istart = iend;
            iend = p[c + 1];

            for (size_t e = istart; e < iend; ++e) {
                r = i[e];
                out[c + r * ncol] = x[e];
            }
        }
}
    }
}



template <typename XITER, typename IITER, typename PITER, 
    typename IT2,
    typename XITER2, typename IITER2, typename PITER2>
extern IT2 csc_rbind(
    std::vector<XITER> const & xvecs, 
    std::vector<IITER> const & ivecs, 
    std::vector<PITER> const & pvecs, 
    std::vector<IT2> const & nrows, 
    std::vector<IT2> const & ncols, 
    XITER2 ox,
    IITER2 oi,
    PITER2 op,
    int const & threads) {

    using IT = typename std::iterator_traits<IITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using iPT = typename std::conditional<std::is_same<PT, double>::value, size_t, int>::type;

    // Rprintf("1\n");
    int n_vecs = nrows.size();

    // check to see if all inputs have the same number of columns.
    IT ncol = ncols[0];
    for (int i = 1; i < n_vecs; ++i) {
        if (ncol != ncols[i]) return 0; // bad input.
    }

    // compute p offsets, 1 per col.
    std::vector<iPT> p_offsets((ncol + 1), 0);
    // Rprintf("2.  ncol = %d, p0 %d, p1 %d, p2799 %d, p2800 %d\n", ncol, p_offsets[0], p_offsets[1], p_offsets[2799], p_offsets[2800]);

    // count the number of non-zero elements per column
    PITER piter = pvecs[0];
    for (IT i = 0; i <= ncol; ++i) {
        iPT offset = 0;
        for (int j = 0; j < n_vecs; ++j) {
            piter = pvecs[j]; piter += i;
            // if ((i < 10) || ( i > 2790)) Rprintf("3. col %d, vec %d, offset %d\n", i, j, static_cast<iPT>(*piter));
            offset += static_cast<iPT>(*piter);
        }
        p_offsets[i] = offset;
        // if ((i < 10) || ( i > 2790)) Rprintf("3. col %d  offset %d, p %d \n", i, offset, p_offsets[i]);
    }
    // Rprintf("3. cols %d  p0 %d, p1 %d, p 2799 %d p2800 %d\n", ncol, p_offsets[0],  p_offsets[1],  p_offsets[2799],  p_offsets[2800]);

    // compute row offsets, 1 per vec, to added to i.
    std::vector<IT> i_offsets(n_vecs+1);  // as many as there are ivecs + 1.
    i_offsets[0] = 0;
    for (int i = 0; i < n_vecs; ++i) {
        i_offsets[i+1] = i_offsets[i] + nrows[i];
    }
    IT nr = i_offsets[n_vecs];
    // size_t nelem  = p_offsets[ncol];

    // using temporary data structure because of omp and R stack imbalance issue.  okay here since we have lots of random writes so rather have the multithreading benefits.
    // Rprintf("4.  nr %d, nz %d\n", nr, nelem);

    // copy data over to the aggregated output
    // iPT nz = 0;
    // auto oxptr = ox;
    // auto oiptr = oi;
    for (int j = 0; j < n_vecs; ++j) {

        XITER ix = xvecs[j];
        IITER ii = ivecs[j];
        PITER ip = pvecs[j];
        // Rprintf("vector %d\n", j);
        for (IT i = 0; i < ncol; ++i) {
            iPT in_start = static_cast<iPT>(*ip);
            ++ip;
            iPT in_end = static_cast<iPT>(*ip);
            size_t out_start = p_offsets[i];

            // Rprintf("vector %d, col %d old poffset %d \n", j, i, out_start);

            // iterate over none zeros. 
            // ox is iterator, not const_iterator.  oxptr does not have an assignment operator
            // iterator has a const member so assignment operator does not work.
            // alternative is to create one.
            XITER2 oxptr(ox + out_start);
            IITER2 oiptr(oi + out_start);
            // xvecs is const & vector.  xvecs[j] is const const_iterator.  const_iterator::operator+ is not const,
            // so xvecs[j] + in_start  does not work.  assignment operator seems to work, so assign first, then using += works.  
            ix = xvecs[j]; ix += in_start;   
            ii = ivecs[j]; ii += in_start;
            for (; in_start < in_end; ++in_start, ++out_start, ++ix, ++ii, ++oxptr, ++oiptr) {
                *oxptr = *ix;  // copy the input x
                *oiptr = *ii + i_offsets[j];
            }
            // Rprintf("vector %d, col %d new poffset %d \n", j, i, out_start);

            p_offsets[i] = out_start;  // set the new position.
        }
    }
    // copy p_offset
    *op = 0;
    // Rprintf("5\n");

#pragma omp parallel num_threads(threads)
{   
    int tid = omp_get_thread_num();
    size_t block = ncol / threads;   // each thread handles a block of rows.
    size_t rem = ncol - threads * block;
    size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

    PITER2 opptr = op + offset + 1;
    for (size_t i = offset; i < end; ++i, ++opptr) {
        // iterator's operator* returns a proxy.  This is the value_type of r_vector<>'s iterator. (const_iterator has T as datatype)
        // creating an array of proxy is not appropriate.  hardcode the type as size_t and let c++ cast automatically.
        (*opptr) = p_offsets[i];
    }
}
    // Rprintf("6\n");

    return nr;
}



template <template <typename...> class VEC, typename XVEC, typename IVEC, typename PVEC, 
    typename COUNTVEC,
    typename XVEC2, typename IVEC2, typename PVEC2>
extern int csc_rbind_vec(
    VEC<XVEC> const & xvecs, 
    VEC<IVEC> const & ivecs, 
    VEC<PVEC> const & pvecs, 
    COUNTVEC const & nrows, 
    COUNTVEC const & ncols, 
    XVEC2 & ox,
    IVEC2 & oi,
    PVEC2 & op,
    int const & threads) {

    int n_vecs = nrows.size();

    // check to see if all inputs have the same number of columns.
    size_t ncol = ncols[0];
    for (int i = 1; i < n_vecs; ++i) {
        if (ncol != static_cast<size_t>(ncols[i])) return 0; // bad input.
    }

    // compute p offsets, 1 per col.
    std::vector<size_t> p_offsets(ncol + 1, 0);

    // global offset = sum of per vec offsets
    // a, a+b, a+b+c
    // d, d+e, d+e+f
    // a+d,  a+d+b+e, ...
    // compute per-column p offset in output.
#pragma omp parallel num_threads(threads)
{   
    int tid = omp_get_thread_num();
    size_t block = (ncol+1) / threads;   // each thread handles a block of rows.
    size_t rem = (ncol+1) - threads * block;
    size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

    size_t off;
    for (size_t i = offset; i < end; ++i) {
        off = 0;
        for (int j = 0; j < n_vecs; ++j) {
            off += static_cast<size_t>(pvecs[j][i]);
        }
        p_offsets[i] = off;
    }
}

    // compute row offsets, 1 per vec, to added to i.
    std::vector<size_t> i_offsets(n_vecs+1);  // as many as there are ivecs + 1.
    i_offsets[0] = 0;
    for (int i = 0; i < n_vecs; ++i) {
        i_offsets[i+1] = i_offsets[i] + nrows[i];
        // Rprintf("i (row) offset for matrix %d = %d\n", i, i_offsets[i]);
    }
    size_t nr = i_offsets[n_vecs];
    // size_t nelem  = p_offsets[ncol];

     // copy data over to the aggregated output
    // size_t nz = 0;
    for (int j = 0; j < n_vecs; ++j) {

        XVEC const & ix = xvecs[j];
        IVEC const & ii = ivecs[j];
        PVEC const & ip = pvecs[j];
        auto off = i_offsets[j];

        for (size_t i = 0; i < ncol; ++i) {
            size_t in_start = static_cast<size_t>(ip[i]);
            size_t in_end = static_cast<size_t>(ip[i+1]);
            size_t out_start = p_offsets[i];

            // iterate over none zeros. 
            // ox is iterator, not const_iterator.  oxptr does not have an assignment operator
            // iterator has a const member so assignment operator does not work.
            // alternative is to create one.
            // xvecs is const & vector.  xvecs[j] is const const_iterator.  const_iterator::operator+ is not const,
            // so xvecs[j] + in_start  does not work.  assignment operator seems to work, so assign first, then using += works.  
            for (; in_start < in_end; ++in_start, ++out_start) {
                ox[out_start] = ix[in_start];  // copy the input x
                oi[out_start] = ii[in_start] + off;
            }

            p_offsets[i] = out_start;  // set the new p offsets after each vector
        }
    }
    // copy p_offset
    op[0] = 0;
// #pragma omp parallel num_threads(threads)
// {   
//     int tid = omp_get_thread_num();
//     size_t block = ncol / threads;   // each thread handles a block of rows.
//     size_t rem = ncol - threads * block;
//     size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
//     int nid = tid + 1;
//     size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);
    
//     for (size_t i = offset; i < end; ++i) {
//         // iterator's operator* returns a proxy.  This is the value_type of r_vector<>'s iterator. (const_iterator has T as datatype)
//         // creating an array of proxy is not appropriate.  hardcode the type as size_t and let c++ cast automatically.
//         op[i + 1] = p_offsets[i];
//     }
// }
    op[0] = 0;
    for (size_t i = 0; i < ncol; ++i) {
        // iterator's operator* returns a proxy.  This is the value_type of r_vector<>'s iterator. (const_iterator has T as datatype)
        // creating an array of proxy is not appropriate.  hardcode the type as size_t and let c++ cast automatically.
        op[i + 1] = p_offsets[i];
    }

    return nr;
}

template <typename XITER, typename IITER, typename PITER, typename IT2,
    typename XITER2, typename IITER2, typename PITER2>
extern IT2 csc_cbind(
    std::vector<XITER> const & xvecs, 
    std::vector<IITER> const & ivecs, 
    std::vector<PITER> const & pvecs, 
    std::vector<IT2> const & nrows, 
    std::vector<IT2> const & ncols, 
    XITER2 ox,
    IITER2 oi,
    PITER2 op,
    int const & threads) {

    using IT = typename std::iterator_traits<IITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using iPT = typename std::conditional<std::is_same<PT, double>::value, size_t, int>::type;

    // data in CSC format. row elements are consecutive, columns are consecutive blocks.
    // rbind requires reconstructing the elements.

    int n_vecs = nrows.size();

    // check to see if all inputs have the same number of rows.
    IT nrow = nrows[0];
    for (int i = 1; i < n_vecs; ++i) {
        if (nrow != nrows[i]) return 0; // bad input.
    }

    // compute p offsets, 1 per vec.
    std::vector<size_t> p_offsets(n_vecs+1);
    std::vector<IT>  c_offsets(n_vecs+1);
    p_offsets[0] = 0;
    c_offsets[0] = 0;
    IT nc = 0;
    for (int i = 0; i < n_vecs; ++i) {
        auto ip = pvecs[i];
        nc = ncols[i];
        c_offsets[i+1] = c_offsets[i] + nc;
        p_offsets[i+1] = p_offsets[i] + static_cast<size_t>(*(ip + nc));
    }
    IT ncol = c_offsets[n_vecs];
    // size_t nelem = p_offsets[n_vecs];

    // cbind would copy big blocks of data, so more memory coherence.
    // parallelize by vecs.

    // int nt = std::min(threads, n_vecs);

    // copy data over to the aggregated output
    // omp causes stack imbalance.   
    iPT nz = 0;
    size_t poff;
    IT coff = 0;
    // auto oxptr = ox;
    // auto oiptr = oi;
    // auto opptr = op;
    for (int j = 0; j < n_vecs; ++j) {

        XITER ix = xvecs[j];
        IITER ii = ivecs[j];
        PITER ip = pvecs[j];

        // copy the x and i vectors
        nc = ncols[j];
        nz = *(ip + nc);
        poff = p_offsets[j];
        XITER2 oxptr(ox + poff);
        IITER2 oiptr(oi + poff);
        
        for (iPT i = 0; i < nz; ++i, ++oxptr, ++oiptr, ++ix, ++ii) {
            *oxptr = *ix;
            *oiptr = *ii;
        }

        coff = c_offsets[j];
        PITER2 opptr(op + coff);
        // copy the p vector
        for (IT i = 0; i < nc; ++i, ++ip, ++opptr) {
            (*opptr) = static_cast<PT>(static_cast<size_t>(*ip) + poff);
        }
    }
    *(op + ncol) = p_offsets[n_vecs];

    return ncol;
}

// r_vector would do a shallow duplicate.   raw pointer would be a copy of pointer only.  
//  c++ vector may involve a copy, but for vector of vector, operator[] returns a reference, 
//  so wthin this funciton we should get a ref to a vector for xvecs[], ivecs[], and pvecs[]
//  
template <template <typename...> class VEC, typename XVEC, typename IVEC, typename PVEC, typename COUNTVEC,
    typename XVEC2, typename IVEC2, typename PVEC2>
extern int csc_cbind_vec(
    VEC<XVEC> const & xvecs, 
    VEC<IVEC> const & ivecs, 
    VEC<PVEC> const & pvecs, 
    COUNTVEC const & nrows, 
    COUNTVEC const & ncols, 
    XVEC2 & ox,
    IVEC2 & oi,
    PVEC2 & op,
    int const & threads) {

    // data in CSC format. row elements are consecutive, columns are consecutive blocks.
    // rbind requires reconstructing the elements.

    int n_vecs = xvecs.size();

    // check to see if all inputs have the same number of rows.
    size_t nrow = nrows[0];
    for (int i = 1; i < n_vecs; ++i) {
        if (nrow != static_cast<size_t>(nrows[i])) return 0; // bad input.
    }

    // compute pointer and column offsets.  as prefix sum.
    std::vector<size_t> p_offsets(n_vecs+1, 0);  // pointer offsets
    std::vector<size_t> c_offsets(n_vecs+1, 0);   // column offsets.  has to be same type as IT2
    size_t nc = 0;
    for (int i = 0; i < n_vecs; ++i) {
        nc = ncols[i];
        c_offsets[i+1] = c_offsets[i] + nc;
        p_offsets[i+1] = p_offsets[i];
        p_offsets[i+1] += static_cast<size_t>(pvecs[i][nc]);
        // Rprintf("%d of %d : NC: %d  coff %d, poff %d\n", i, n_vecs, nc, c_offsets[i], p_offsets[i]);
    }
    size_t ncol = c_offsets[n_vecs];   // total column count
    size_t nelem = p_offsets[n_vecs];   // total element count

    // cbind would copy big blocks of data, so more memory coherence.
    // parallelize by vecs.

    // int nt = std::min(threads, n_vecs);

    // copy data over to the aggregated output
    // omp causes stack imbalance.   
    size_t nz = 0;
    size_t poff = 0;
    size_t coff = 0;

    for (int j = 0; j < n_vecs; ++j) {

        XVEC const & ix = xvecs[j];
        IVEC const & ii = ivecs[j];
        PVEC const & ip = pvecs[j];

        nc = ncols[j];
        nz = ip[nc];  // get the currnet vec's non-zero count

        // copy the p vector        
        coff = c_offsets[j];
        poff = p_offsets[j];
        // Rprintf("%d of %d : coff %d, poff %d\n", j, n_vecs, coff, poff);

        // copy the p vector
        for (size_t i = 0; i < nc; ++i, ++coff) {
            op[coff] = ip[i] + poff;
            // Rprintf("\t%d ip: %d, poff: %d, coff %d result %d\n", i, ip[i], poff, coff, offset);
            // op[coff] = offset;
        }

        // copy the x and i vectors
        for (size_t i = 0; i < nz; ++i, ++poff) {
            ox[poff] = ix[i];
            oi[poff] = ii[i];
        }

    }
    op[ncol] = nelem;
    
    // Rprintf("nelem: %d, ncol: %d\n", nelem, ncol);

    return ncol;
}

// csc
template <typename OITER, typename XITER, typename PITER, typename IT2>
extern void csc_colsums_iter(
    XITER x, 
    PITER p, 
    IT2 const & ncol, 
    OITER out,
    int const & threads) {

    using XT = typename std::iterator_traits<XITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using iPT = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

    if (threads == 1) {
        auto pend = p + ncol + 1;
        auto xptr = x;
        auto xend = x;

        iPT start2 = *p;
        ++p;
        iPT end2;

        auto optr = out;
        for (size_t c = 0; c < ncol; ++c, ++optr) {*optr = 0; }
        // std::fill(out, out + ncol, 0);   // absolutely needed here.
        for (; p != pend; ++p, ++out) {
            end2 = *p;

            XT sum = 0;
            xptr = x + start2;
            xend = x + end2;
            for (; xptr != xend; ++xptr) {
                sum += *xptr;
            }

            *out = sum;
            start2 = end2;
        }

    } else {

#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        auto pptr = p + offset;
        auto pend = p + end;
        OITER optr(out + offset);
        auto xptr = x;
        auto xend = x + *pptr;

        auto oend = out + end;
        for (; optr != oend; ++optr) { *optr = 0; }
        // std::fill(out + offset, out + end, 0);
        optr = out; optr += offset;
        for (; pptr != pend; ++optr) {
            xptr = xend;
            ++pptr;
            xend = x + *pptr;

            XT sum = 0;
            for (; xptr != xend; ++xptr) {
                sum += *xptr;
            }

            *optr = sum;
        }
}
    }
}


// csc
template <typename OVEC, typename XVEC, typename PVEC, typename IT2>
extern void csc_colsums_vec(
    XVEC const & x, 
    PVEC const & p, 
    IT2 const & ncol,
    OVEC & out,   // should be sized already.
    int const & threads) {

    // rely on operator[]
    if (threads == 1) {
        size_t xend = p[0];  // first entry will be star 
        size_t xstart = xend; 
        auto sum = x[0];   // assign just to get hte type

        for (IT2 c = 0; c < ncol; ++c) {
            xstart = xend;
            xend = p[c+1];

            sum = 0;
            for (; xstart < xend; ++xstart) {
                sum += x[xstart];
            }

            out[c] = sum;
        }

    } else {

#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);


        size_t xend = p[offset];  // first entry will be star 
        size_t xstart = xend; 
        auto sum = x[offset];   // assign just to get hte type

        for (IT2 c = offset; c < end; ++c) {
            xstart = xend;
            xend = p[c+1];

            sum = 0;
            for (; xstart < xend; ++xstart) {
                sum += x[xstart];
            }

            out[c] = sum;
        }
}
    }
}

// csc
template <typename OITER, typename XITER, typename IITER, typename IT2, typename IT>
extern void csc_rowsums_iter(
    XITER x, 
    IITER i, 
    IT2 const & nrow, IT const & nzcount, 
    OITER out,
    int const & threads) {

    using XT = typename std::iterator_traits<XITER>::value_type;

    if (threads == 1) {
        // std::vector<XT> sums(nrow, 0);

        size_t r;
        OITER oend = out + nrow;
        for (OITER optr = out; optr != oend; ++optr) {
            *optr = 0;
        }
        IITER end = i + nzcount;
        for (; i != end; ++i, ++x) {
            r = *i;
            *(out + r) += *x;
        }

        // std::copy(sums.begin(), sums.end(), out);
    } else {

        // alloc temp storage.
        std::vector<std::vector<XT>> sums(threads);

#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nzcount / threads;
        size_t rem = nzcount - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        sums[tid].resize(nrow, 0);
        std::vector<XT> & s = sums[tid];

        // IT2 r;
        auto iptr = i + offset;
        auto iend = i + end;
        auto xptr = x + offset;
        for (; iptr != iend; ++iptr, ++xptr) {
            s[*iptr] += *xptr;
        }
}


#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nrow / threads;
        size_t rem = nrow - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        OITER optr(out + offset);
        for (; offset < end; ++offset, ++optr) {
            XT sum = 0;
            for (int t = 0; t < threads; ++t) {
                sum += sums[t][offset];
            }

            *optr = sum;
        }
}

    }

}



// csc
template <typename OVEC, typename XVEC, typename IVEC, typename IT2>
extern void csc_rowsums_vec(
    XVEC const & x, 
    IVEC const & i, 
    IT2 const & nrow,
    OVEC & out,
    int const & threads) {

    size_t nzcount = x.size();

    if (threads == 1) {
        // have to initialize.
        for (IT2 r = 0; r < nrow; ++r) {
            out[r] = 0;
        }
        for (size_t e = 0; e < nzcount; ++e) {
            out[i[e]] += x[e];
        }

    } else {
        using XT = typename std::remove_cv<typename std::remove_reference<decltype(x[0])>::type>::type;
        // alloc temp storage.
        std::vector<std::vector<XT>> sums(threads);


#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nzcount / threads;
        size_t rem = nzcount - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        sums[tid].resize(nrow, 0); // initialize
        std::vector<XT> & s = sums[tid];

        for (size_t e = offset; e < end; ++e) {
            s[i[e]] += x[e];
        }
}


#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = nrow / threads;
        size_t rem = nrow - threads * block;
        size_t offset = tid * block + (static_cast<size_t>(tid) > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (static_cast<size_t>(nid) > rem ? rem : nid);

        for (size_t r = offset; r < end; ++r) {
            XT sum = 0;
            for (int t = 0; t < threads; ++t) {
                sum += sums[t][r];
            }

            out[r] = sum;
        }
}

    }

}
