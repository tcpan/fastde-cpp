#include "utils_sparsemat.hpp"

#include "cpp11/doubles.hpp"
#include "cpp11/integers.hpp"

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


// #include <atomic>
// #include <parallel/algorithm>
// VERY SLOW
// // NOTe:  there is no formal definition of sparse matrix.
// // input is column major, so i has the row ids, and p is per column.
// template <typename XT, typename IT, typename PT, typename IT2>
// extern cpp11::writable::list _sp_transpose_sort(
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
//     size_t offset = tid * block + (tid > rem ? rem : tid);
//     int nid = tid + 1;
//     size_t end = nid * block + (nid > rem ? rem : nid);

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
//     size_t offset = tid * block + (tid > rem ? rem : tid);

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
//     size_t offset = tid * block + (tid > rem ? rem : tid);
//     int nid = tid + 1;
//     size_t end = nid * block + (nid > rem ? rem : nid);

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



// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template <typename XT, typename IT, typename PT, typename IT2>
extern cpp11::writable::list _sp_transpose_par(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, int const & threads) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

  std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;


  start = std::chrono::steady_clock::now();


    // either:  per thread summary,   this would still use less memory than sortiing the whole thing.
    // bin by row in random (thread) order, then sort per row -  this would be n log n - n log t
    
    PT2 nelem = x.size();

    // empty output 
    cpp11::writable::r_vector<XT> tx(nelem); 
    cpp11::writable::r_vector<IT> ti(nelem);   // as many as there are values 
    cpp11::writable::r_vector<PT> tp(nrow + 1);       // number of rows + 1.

  Rprintf("[TIME] sp_transpose_par 1 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

    // set up per thread offsets.
    // should have
    //      r0      r1      r2      r3  ...
    // t0   0       s1      s3
    // t1   c1      s1+d1
    // t2   c1+c2   s1+d1+d2
    // ...
    // tn   s1      s2
    std::vector<std::vector<PT>> lps(threads + 1);


#pragma omp parallel num_threads(threads)
{   
    int tid = omp_get_thread_num();
    lps[tid].resize(nrow + 1);
    std::fill(lps[tid].begin(), lps[tid].end(), 0);
}
    lps[threads].resize(nrow + 1);
    std::fill(lps[threads].begin(), lps[threads].end(), 0);

  Rprintf("[TIME] sp_transpose_par 2 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

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
    size_t offset = tid * block + (tid > rem ? rem : tid);
    int nid = tid + 1;
    size_t end = nid * block + (nid > rem ? rem : nid);

    for (; offset != end; ++offset) {
        ++lps[tid+1][static_cast<IT2>(i[offset])];
    }
}
  Rprintf("[TIME] sp_transpose_par 3 Elapsed(ms)= %f\n", since(start).count());

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

    for (int t = 1; t <= threads; ++t) {  // linear scan, for hardware prefetching.
        for (size_t r = offset; r != end; ++r) {
            lps[t][r] += lps[t-1][r];
        }
    }
    // at the end, lps[thread] has total counts per row.
}
      Rprintf("[TIME] sp_transpose_par 4 Elapsed(ms)= %f\n", since(start).count());

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
        tp[r+1] = lps[0][r+1];
    }
    tp[0] = lps[0][0];
  Rprintf("[TIME] sp_transpose_par 5 Elapsed(ms)= %f\n", since(start).count());

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
  Rprintf("[TIME] sp_transpose_par 6 Elapsed(ms)= %f\n", since(start).count());

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
    auto pptr = std::upper_bound(p.begin(), p.end(), offset);
    IT2 cid = std::distance(p.begin(), pptr) - 1;
    PT2 pos;
    
    for (; offset < end; ++offset) {
        rid = i[offset];   // current row id (starts with 0)
        val = x[offset];   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; offset >= static_cast<PT2>(p[cid+1]); ++cid);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = lps[tid][rid];  // where to insert the data
        tx[pos] = val;  // place the data
        ti[pos] = cid;  // place the row id (original col id. 0-based)
        ++lps[tid][rid];  // update the offset - 1 space consumed.
    }
}
  Rprintf("[TIME] sp_transpose_par 7 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

    // ======= return
    cpp11::named_arg _tx("x"); _tx = tx;
    cpp11::named_arg _ti("i"); _ti = ti;
    cpp11::named_arg _tp("p"); _tp = tp;
    cpp11::writable::list out( { _tx, _ti, _tp} );

      Rprintf("[TIME] sp_transpose_par 8 Elapsed(ms)= %f\n", since(start).count());

    return out;

}




// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
template <typename XT, typename IT, typename PT, typename IT2>
extern cpp11::writable::list _sp_transpose(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, IT2 const & nrow, IT2 const & ncol, int const & threads) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

    // input
    
    // empty output 
    cpp11::writable::r_vector<XT> tx(x.size()); 
    cpp11::writable::r_vector<IT> ti(x.size());   // as many as there are values 
    cpp11::writable::r_vector<PT> tp(nrow + 1);       // number of rows + 1.
    // init tp
    std::fill(tp.begin(), tp.end(), 0);
        
    // ======= do the transpose.
    
    // do the swap.  do random memory access instead of sorting.
    // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
    // 2. compute exclusive prefix sum as offsets in x.
    // 3. use p to get range of elements in x belonging to the same column, scatter to new x
    //     and increment offset in p.
    // 4. shift p to the right by 1 element, and set p[0] = 0
    // step 1
    auto iptr = i.begin();
    auto i_end = i.end();  // remember i is 0-based.
    for (; iptr != i_end; ++iptr) {
        ++tp[static_cast<IT2>(*iptr) + 1];
    }
    // step 2 - create max offset + 1 for each transposed row ( == new column)
    auto tncol = nrow;
    IT2 cid = 1;
    for (; cid <= tncol; ++cid) {
        tp[cid] += tp[cid - 1];
    }
    // step 3  scatter
    IT2 rid;   // column id needs to start with 0.  row ids start with 0
    auto pptr = p.begin() + 1;  // compare to end of col ptr
    cid = 0;
    XT val;
    PT2 nelem = x.size();
    PT2 pos;
    decltype(nelem) e = 0;
    
    for (PT2 e = 0; e < nelem; ++e) {
        rid = i[e];   // current row id (starts with 0)
        val = x[e];   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; e >= static_cast<PT2>(*pptr); ++cid, ++pptr);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = tp[rid];  // where to insert the data
        tx[pos] = val;  // place the data
        ti[pos] = cid;  // place the row id (original col id. 0-based)
        ++tp[rid];  // update the offset - 1 space consumed.
    }
    // step 4
    PT temp = tp[0];
    PT temp2;
    for (IT2 cid = 1; cid <= tncol; ++cid) {
        temp2 = tp[cid];
        tp[cid] = temp;
        temp = temp2;
    }
    tp[0] = 0;

    // ======= return
    cpp11::named_arg _tx("x"); _tx = tx;
    cpp11::named_arg _ti("i"); _ti = ti;
    cpp11::named_arg _tp("p"); _tp = tp;
    return cpp11::writable::list( { _tx, _ti, _tp} );
}


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XT, typename IT, typename PT, typename IT2, typename PT2>
extern void _sp_transpose_par(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol,
    std::vector<XT> & tx, 
    std::vector<IT2> & ti, 
    std::vector<PT2> & tp, 
    int const & threads) {

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
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;

  start = std::chrono::steady_clock::now();

    // empty output 
    tx.clear();    tx.resize(x.size()); 
    ti.clear();    ti.resize(x.size());   // as many as there are values 
    tp.clear();    tp.resize(nrow + 1);       // number of rows + 1.
  Rprintf("[TIME] sp_transpose_par alloc output Elapsed(ms)= %f\n", since(start).count());

    _sp_transpose_par(x, i, p, nrow, ncol, tx.data(), ti.data(), tp.data(), threads);    


}


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XT, typename IT, typename PT, typename IT2, typename PT2>
extern void _sp_transpose_par(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol,
    XT * tx, 
    IT2 * ti, 
    PT2 * tp, 
    int const & threads) {

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
    std::chrono::time_point<std::chrono::steady_clock, std::chrono::duration<double>> start;

  start = std::chrono::steady_clock::now();

    PT2 nelem = x.size();

    // assume all allocated properly

  Rprintf("[TIME] sp_transpose_par 1 Elapsed(ms)= %f\n", since(start).count());

  start = std::chrono::steady_clock::now();

    // set up per thread offsets.
    // should have
    //      r0      r1      r2      r3  ...
    // t0   0       s1      s3
    // t1   c1      s1+d1
    // t2   c1+c2   s1+d1+d2
    // ...
    // tn   s1      s2
    std::vector<std::vector<PT>> lps(threads + 1);
#pragma omp parallel num_threads(threads)
{   
    int tid = omp_get_thread_num();
    lps[tid].resize(nrow + 1);
    std::fill(lps[tid].begin(), lps[tid].end(), 0);
}
    lps[threads].resize(nrow + 1);
    std::fill(lps[threads].begin(), lps[threads].end(), 0);

    // ======= do the transpose.
      Rprintf("[TIME] sp_transpose_par 2 Elapsed(ms)= %f\n", since(start).count());

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

    for (; offset != end; ++offset) {
        ++lps[tid+1][static_cast<IT2>(i[offset])];
    }
}
  Rprintf("[TIME] sp_transpose_par 3 Elapsed(ms)= %f\n", since(start).count());

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

    for (int t = 1; t <= threads; ++t) {  // linear scan, for hardware prefetching.
        for (size_t r = offset; r != end; ++r) {
            lps[t][r] += lps[t-1][r];
        }
    }
    // at the end, lps[thread] has total counts per row.
}
      Rprintf("[TIME] sp_transpose_par 4 Elapsed(ms)= %f\n", since(start).count());

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
        tp[r+1] = lps[0][r+1];
    }
    tp[0] = lps[0][0];

  Rprintf("[TIME] sp_transpose_par 5 Elapsed(ms)= %f\n", since(start).count());

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
  Rprintf("[TIME] sp_transpose_par 6 Elapsed(ms)= %f\n", since(start).count());

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
    auto pptr = std::upper_bound(p.begin(), p.end(), offset);
    IT2 cid = std::distance(p.begin(), pptr) - 1;
    PT2 pos;
    
    for (; offset < end; ++offset) {
        rid = i[offset];   // current row id (starts with 0)
        val = x[offset];   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; offset >= static_cast<PT2>(p[cid+1]); ++cid);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = lps[tid][rid];  // where to insert the data
        tx[pos] = val;  // place the data
        ti[pos] = cid;  // place the row id (original col id. 0-based)
        ++lps[tid][rid];  // update the offset - 1 space consumed.
    }
}
  Rprintf("[TIME] sp_transpose_par 7 Elapsed(ms)= %f\n", since(start).count());


}




// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XT, typename IT, typename PT, typename IT2, typename PT2>
extern void _sp_transpose(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    std::vector<XT> & tx, 
    std::vector<IT2> & ti, 
    std::vector<PT2> & tp, 
    int const & threads) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // input

    // empty output 
    tx.clear(); tx.resize(x.size()); 
    ti.clear(); ti.resize(x.size());   // as many as there are values 
    tp.clear(); tp.resize(nrow + 1);       // number of rows + 1.

    _sp_transpose(x, i, p, nrow, ncol, tx.data(), ti.data(), tp.data(), threads);

}


// NOTe:  there is no formal definition of sparse matrix.
// input is column major, so i has the row ids, and p is per column.
// direct to stl vectors
template <typename XT, typename IT, typename PT, typename IT2, typename PT2>
extern void _sp_transpose(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, 
    XT * tx, 
    IT2 * ti, 
    PT2 * tp, 
    int const & threads) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    // input

    // assume output pointers are all allocated
    // init tp
    memset(tp, 0, (nrow + 1) * sizeof(PT2));
        
    // ======= do the transpose.
    
    // do the swap.  do random memory access instead of sorting.
    // 1. iterate over i to get row (tcol) counts, store in new p[1..nrow].   these are offsets in new x
    // 2. compute exclusive prefix sum as offsets in x.
    // 3. use p to get range of elements in x belonging to the same column, scatter to new x
    //     and increment offset in p.
    // 4. shift p to the right by 1 element, and set p[0] = 0
    // step 1
    auto iptr = i.begin();
    auto i_end = i.end();  // remember i is 0-based.
    for (; iptr != i_end; ++iptr) {
        ++tp[static_cast<IT2>(*iptr) + 1];
    }
    // step 2 - create max offset + 1 for each transposed row ( == new column)
    auto tncol = nrow;
    IT2 cid = 1;
    for (; cid <= tncol; ++cid) {
        tp[cid] += tp[cid - 1];
    }
    // step 3  scatter
    IT2 rid;   // column id needs to start with 0.  row ids start with 0
    auto pptr = p.begin() + 1;  // compare to end of col ptr
    cid = 0;
    XT val;
    PT2 nelem = x.size();
    PT2 pos;
    decltype(nelem) e = 0;
    
    for (PT2 e = 0; e < nelem; ++e) {
        rid = i[e];   // current row id (starts with 0)
        val = x[e];   // current value
        // if the current element pos reaches first elem of next column (*pptr),
        // then go to next column (increment cid and pptr).
        for (; e >= static_cast<PT2>(*pptr); ++cid, ++pptr);  // current column id

        // now copy and update.
        // curr pos is the offset for the transposed row (new col), in tp.
        // note we are using tp array to track current offset.
        pos = tp[rid];  // where to insert the data
        tx[pos] = val;  // place the data
        ti[pos] = cid;  // place the row id (original col id. 0-based)
        ++tp[rid];  // update the offset - 1 space consumed.
    }
    // step 4
    PT temp = tp[0];
    PT temp2;
    for (IT2 cid = 1; cid <= tncol; ++cid) {
        temp2 = tp[cid];
        tp[cid] = temp;
        temp = temp2;
    }
    tp[0] = 0;

}


// no names.
template <typename OUT, typename XT, typename IT, typename PT, typename IT2>
extern OUT _sp_to_dense(cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, cpp11::r_vector<PT> const & p, 
    IT2 const & nrow, IT2 const & ncol, int const & threads) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
    
    using PT2 = typename std::conditional<std::is_same<PT, double>::value, long, int>::type;

    // ======= create new output and initialize
    OUT out(nrow, ncol);
    auto slice_end = out.end();
    for (auto slice = out.begin(); slice != slice_end; ++slice) {
        std::fill((*slice).begin(), (*slice).end(), 0);
    }

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    if (threads == 1) {
        IT2 r;
        PT2 istart, iend;
        for (size_t c = 0; c < ncol; ++c) {
            istart = p[c];
            iend = p[c+1];

            for (; istart < iend; ++istart) {
                r = i[istart];
                out(r, c) = x[istart];
            }
        }

    } else {
    // now iterate and fill   
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
            istart = p[offset];
            iend = p[offset+1];

            for (; istart < iend; ++istart) {
                r = i[istart];
                out(r, offset) = x[istart];
            }
        }
}

    }
    // now iterate and fill   

    return out;
}


template <typename OUT, typename XT, typename IT, typename PT, typename IT2>
extern OUT _sp_to_dense_transposed(
    cpp11::r_vector<XT> const & x, 
    cpp11::r_vector<IT> const & i, 
    cpp11::r_vector<PT> const & p, IT2 const & nrow, IT2 const & ncol, int const & threads) {

    // https://www.r-bloggers.com/2020/03/what-is-a-dgcmatrix-object-made-of-sparse-matrix-format-in-r/
    // ======= decompose the input matrix in CSC format, S4 object with slots:
    // i :  int, row numbers, 0-based.
    // p :  int, p[i] is the position offset in x for row i.  i has range [0-r] inclusive.
    // x :  numeric, values
    // Dim:  int, 2D, sizes of full matrix
    // Dimnames:  2D, names.
    // factors:  ignore.
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
    OUT out(ncol, nrow);
    auto slice_end = out.end();
    for (auto slice = out.begin(); slice != slice_end; ++slice) {
        std::fill((*slice).begin(), (*slice).end(), 0);
    }

    // Rprintf("Sparse DIM: samples %lu x features %lu, non-zeros %lu\n", in.get_ncol(), in.get_nrow(), in.get_nelem()); 

    if (threads == 1) {
    // now iterate and fill   
        IT r;
        PT2 istart, iend;
        for (size_t c = 0; c < ncol; ++c) {  // iterate by source.
            istart = p[c];
            iend = p[c+1];

            for (; istart < iend; ++istart) {
                r = i[istart];
                out(c, r) = x[istart];
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
            istart = p[offset];
            iend = p[offset+1];

            for (; istart < iend; ++istart) {
                r = i[istart];
                out(offset, r) = x[istart];
            }
        }
}

    }
 
    return out;
}
