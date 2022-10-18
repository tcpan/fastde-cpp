#include "fastde/sparsemat.hpp"

/*
 * wrapper for R dgCMatrix
 *
 */

#include <vector>
#include <algorithm>
#include <cstring>

#include <omp.h>

#include "utils_data.hpp"
#include "fastde/benchmark_utils.hpp"



// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XITER, typename IITER, typename PITER, typename XITER2, typename IITER2, typename PITER2>
extern void _sp_transpose_par(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol,
    XITER2 tx, 
    IITER2 ti, 
    PITER2 tp, 
    int const & threads) {

    using XT = typename std::iterator_traits<XITER>::value_type;
    using IT = typename std::iterator_traits<IITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    PT2 nelem = *(p + ncol);

    // either:  per thread summary,   this would still use less memory than sortiing the whole thing.
    // bin by row in random (thread) order, then sort per row -  this would be n log n - n log t
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;


  start = std::chrono::steady_clock::now();

    // initialize per thread offsets.
    // should have
    //      r0      r1      r2      r3  ...
    // t0   0       s1      s3
    // t1   c1      s1+d1
    // t2   c1+c2   s1+d1+d2
    // ...
    // tn   s1      s2
    PT * gp = reinterpret_cast<PT *>(malloc((threads + 1) * (nrow + 1)));
    std::vector<PT *> lps(threads + 1);
#pragma omp parallel num_threads(threads)
{   
    int tid = omp_get_thread_num();
    lps[tid] = gp + tid * (nrow + 1);
    memset(lps[tid], 0, (nrow + 1) * sizeof(PT));
}
    lps[threads] = gp + threads * (nrow + 1);
    memset(lps[threads], 0, (nrow + 1) * sizeof(PT));

    // ======= do the transpose.
    // Rprintf("[TIME] sp_transpose_par 2 Elapsed(ms)= %f\n", since(start).count());

    start = std::chrono::steady_clock::now();

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
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    auto istart = i + offset;
    auto iend = i + end;
    for (; istart != iend; ++istart) {
        ++lps[tid+1][static_cast<IT2>(*istart)];
    }
}
//   Rprintf("[TIME] sp_transpose_par 3 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

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
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    for (int t = 0; t < threads; ++t) {  // linear scan, for hardware prefetching.
        for (size_t r = offset; r != end; ++r) {
            lps[t + 1][r] += lps[t][r];
        }
    }
    // at the end, lps[thread] has total counts per row.
}
    //   Rprintf("[TIME] sp_transpose_par 4 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

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
        *(tp + r+1) = lps[0][r+1];
    }
    *tp = lps[0][0];

//   Rprintf("[TIME] sp_transpose_par 5 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

    // step 2.1: add global prefix to local.  do in parallel.  each thread can do independently.
    //      r0      r1      r2      r3  ...
    // t0   0       s1      s3
    // t1   c1      s1+d1
    // t2   c1+c2   s1+d1+d2
    // ...
    // tn   s1      s1+s2
#pragma omp parallel num_threads(threads)
{   
    int tid = omp_get_thread_num();
    for (IT2 r = 0; r < nrow; ++r) {
        lps[tid + 1][r] += lps[0][r];
    }
}
    // per thread we now have the starting offset for writing.
//   Rprintf("[TIME] sp_transpose_par 6 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

    // step 3.  use the per thread offsets to write out.
#pragma omp parallel num_threads(threads)
{   
    int tid = omp_get_thread_num();
    size_t block = nelem / threads;
    size_t rem = nelem - threads * block;
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    IT2 rid;   // column id needs to start with 0.  row ids start with 0
    XT val;
    // need to search for cid based on offset.
    auto pptr = std::upper_bound(p, p + nrow + 1, offset);
    IT2 cid = std::distance(p, pptr) - 1;
    PT2 pos;
    auto iptr = i + offset;
    auto xptr = x + offset;
    
    for (; offset < end; ++offset, ++iptr, ++xptr) {
        rid = *iptr;   // current row id (starts with 0)
        val = *xptr;   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; offset >= static_cast<PT2>(*pptr); ++cid, ++pptr);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = lps[tid][rid];  // where to insert the data
        *(tx + pos) = val;  // place the data
        *(ti + pos) = cid;  // place the row id (original col id. 0-based)
        ++lps[tid][rid];  // update the offset - 1 space consumed.
    }
}
//   Rprintf("[TIME] sp_transpose_par 7 Elapsed(ms)= %f\n", since(start).count());


}

// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XITER, typename IITER, typename PITER, 
    typename XITER2, typename IITER2, typename PITER2>
extern void _sp_transpose(
    XITER x, 
    IITER i, 
    PITER p, 
    IT2 const & nrow, IT2 const & ncol,
    XITER2 tx, 
    IITER2 ti, 
    PITER2 tp, 
    int const & threads) {

    using XT = typename std::iterator_traits<XITER>::value_type;
    using IT = typename std::iterator_traits<IITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // input
    PT2 nelem = *(p + ncol);

    // assume output pointers are all allocated
    // init tp
    std::fill(tp, tp + nrow + 1, 0);
        
    // ======= do the transpose.
    
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
    PT2 pos;
    auto tpptr = tp;
    for (size_t e = 0; e < nelem; ++e, ++i, ++x) {
        rid = *i;   // current row id (starts with 0)
        val = *x;   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; e >= static_cast<PT2>(*pptr); ++cid, ++pptr);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = *(tpptr + rid);  // where to insert the data
        *(tx + pos) = val;  // place the data
        *(ti + pos) = cid;  // place the row id (original col id. 0-based)
        *(tpptr + rid) += 1;  // update the offset - 1 space consumed.
    }
    // step 4
    std::copy_backward(tp, tp + nrow, tp + nrow + 1);
    *tp = 0;
}


// cpp11 r_vector iterator is slow compared to bracket accessor.  
//  It is not clear if this is due to overhead of iterator instantiation or . 
template <typename OITER, typename XITER, typename IITER, typename PITER, typename IT2>
extern void _sp_to_dense(
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
    using XT = typename std::iterator_traits<XITER>::value_type;
    using IT = typename std::iterator_traits<IITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;


    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // ======= create new output and initialize
    std::fill(out, out + (nrow * ncol), 0);

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    if (threads == 1) {
    // now iterate and fill   
        IT r;
        PT2 istart, iend;
        for (size_t c = 0; c < ncol; ++c) {  // iterate by source.
            istart = *(p + c);
            iend = *(p + c + 1);

            for (; istart < iend; ++istart) {
                r = *(i + istart);
                out(r, c) = *(x + istart);
            }
        }
    } else {
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (tid > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (nid > rem ? rem : nid);

        IT2 r;
        PT2 istart, iend;
        for (; offset < end; ++offset) {
            istart = *(p+offset);
            iend = *(p+offset+1);

            for (; istart < iend; ++istart) {
                r = *(i+istart);
                out(r, offset) = *(x+istart);
            }
        }
}

    }
 
}


template <typename OITER, typename XITER, typename IITER, typename PITER, typename IT2>
extern void _sp_to_dense_transposed(
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
    using XT = typename std::iterator_traits<XITER>::value_type;
    using IT = typename std::iterator_traits<IITER>::value_type;
    using PT = typename std::iterator_traits<PITER>::value_type;
    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;


    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // ======= create new output and initialize
    std::fill(out, out + (nrow * ncol), 0);

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    if (threads == 1) {
    // now iterate and fill   
        IT r;
        PT2 istart, iend;
        for (size_t c = 0; c < ncol; ++c) {  // iterate by source.
            istart = *(p + c);
            iend = *(p + c + 1);

            for (; istart < iend; ++istart) {
                r = *(i + istart);
                out(c, r) = *(x + istart);
            }
        }
    } else {
#pragma omp parallel num_threads(threads)
{   
        int tid = omp_get_thread_num();
        size_t block = ncol / threads;
        size_t rem = ncol - threads * block;
        size_t offset = tid * block + (tid > rem ? rem : tid);
        int nid = tid + 1;
        size_t end = nid * block + (nid > rem ? rem : nid);

        IT2 r;
        PT2 istart, iend;
        for (; offset < end; ++offset) {
            istart = *(p+offset);
            iend = *(p+offset+1);

            for (; istart < iend; ++istart) {
                r = *(i+istart);
                out(offset, r) = *(x+istart);
            }
        }
}

    }
 
}
