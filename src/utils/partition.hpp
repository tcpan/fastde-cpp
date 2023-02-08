#pragma once

#include <cstring>
#include <vector>
#include <algorithm>
#include <type_traits>

#ifdef USE_MPI
#include <mpi.h>
#include "splash/utils/mpi_types.hpp"
#endif

namespace utils {

// this file contains 1D and 2D partitioners.
//      start with full partitioners (full 1D, and full 2D rectilinear, fixed size blocks or equipartitioned.)
//      the full partitioner produces a linearized list of partitions
//      when partitioning, the parameter can be specified using a partition object as well.
//      
// subsequent filters can be used to rearrange the linear layout and exclude certain parts.

#define PARTITION_EQUAL 1
#define PARTITION_FIXED 2

#define PARTITION_TILE_DIM 8

// no real need for set operations (union, intersection, difference)
// really only need 1. create across multiple nodes.  2. divide, and 3. get specific one.
// template <typename ST, size_t DIMS>
// class block {
//     public:
//         ST[DIMS] offset;
//         ST[DIMS] size;

//         virtual bool contains(ST[DIMS] pt) const {
//             bool out = true;
//             for (size_t i = 0; i < DIMS; ++i) {
//                 out &= (pt[i] >= offset[i]) && (pt[i] < (offset[i] + size[i])); 
//             }
//             return out;
//         };

//         // overlap?
//         virtual bool overlap(block<ST, DIMS> const & rhs) const {
//             return this->contains(rhs.offset) || rhs.contains(this->offset);
//         }

//         virtual block<ST, DIMS> set_union(block<ST, DIMS> const & rhs) {
//             block<ST, DIMS> out;
//             if (this->overlap(rhs)) {
//                 for ( size_t i = 0; i < DIMS; ++i) {
//                     out.offset[i] = 
//                 }
//             } 
//             return out;
//         };
//         // OR operator. union
//         inline block<ST, DIMS>& operator|=(block<ST, DIMS> const & rhs) { 
//             return this->set_union(rhs); 
//         };
//         // passing lhs by value helps optimize chained a+b+c.  then return results by value (move contructor, elision)
//         friend block<ST, DIMS> operator|(block<ST, DIMS> lhs, block<ST, DIMS> const & rhs) {
//             return lhs.set_union(lhs);
//         }; 

//         virtual block<ST, DIMS>& set_intersection(block<ST, DIMS> const & rhs) = 0;
//         // AND operator.  intersection
//         inline block<ST, DIMS>& operator&=(block<ST, DIMS> const & rhs) {
//             return this->set_intersection(rhs);
//         };
//         // passing lhs by value helps optimize chained a+b+c
//         friend block<ST, DIMS> operator&(block<ST, DIMS> lhs, block<ST, DIMS> const & rhs) {
//             return lhs->set_intersection(rhs);
//         }; 

//         virtual block<ST, DIMS>& set_difference(block<ST, DIMS> const & rhs) = 0;
//         // set difference operator. 
//         inline block<ST, DIMS>& operator-=(block<ST, DIMS> const & rhs) {
//             return this->set_difference(rhs);
//         };
//         // passing lhs by value helps optimize chained a+b+c
//         friend block<ST, DIMS> operator-(block<ST, DIMS> lhs, block<ST, DIMS> const & rhs) {
//             return lhs.set_difference(rhs);
//         }; 

//         // equal partition
//         virtual std::vector<block<ST, DIMS>> set_equal_partition(size_t[DIMS] count) = 0;
//         virtual block<ST, DIMS> set_equal_partition_ith(size_t[DIMS] count, size_t[DIMS] idx) = 0;

//         // fixed partition
//         virtual std::vector<block<ST, DIMS>> set_fixed_partition(ST[DIMS] size) = 0;
//         virtual block<ST, DIMS> set_equal_partition_ith(ST[DIMS] size, size_t[DIMS] idx) = 0;
        
// };


// template<typename ST>
// class block<ST, 1> {
//     public:
//         virtual bool inside(ST[DIMS] pt);

//         // overlap?
//         virtual bool overlap(block<ST, DIMS> const & rhs) = 0;

//         virtual block<ST, DIMS>& set_union(block<ST, DIMS> const & rhs) = 0;
//         // OR operator. union
//         inline block<ST, DIMS>& operator|=(block<ST, DIMS> const & rhs) { 
//             return this->set_union(rhs); 
//         };
//         // passing lhs by value helps optimize chained a+b+c.  then return results by value (move contructor, elision)
//         friend block<ST, DIMS> operator|(block<ST, DIMS> lhs, block<ST, DIMS> const & rhs) {
//             return lhs.set_union(lhs);
//         }; 

//         virtual block<ST, DIMS>& set_intersection(block<ST, DIMS> const & rhs) = 0;
//         // AND operator.  intersection
//         inline block<ST, DIMS>& operator&=(block<ST, DIMS> const & rhs) {
//             return this->set_intersection(rhs);
//         };
//         // passing lhs by value helps optimize chained a+b+c
//         friend block<ST, DIMS> operator&(block<ST, DIMS> lhs, block<ST, DIMS> const & rhs) {
//             return lhs->set_intersection(rhs);
//         }; 

//         virtual block<ST, DIMS>& set_difference(block<ST, DIMS> const & rhs) = 0;
//         // set difference operator. 
//         inline block<ST, DIMS>& operator-=(block<ST, DIMS> const & rhs) {
//             return this->set_difference(rhs);
//         };
//         // passing lhs by value helps optimize chained a+b+c
//         friend block<ST, DIMS> operator-(block<ST, DIMS> lhs, block<ST, DIMS> const & rhs) {
//             return lhs.set_difference(rhs);
//         }; 

//         // equal partition
//         virtual std::vector<block<ST, DIMS>> set_equal_partition(size_t[DIMS] count) = 0;
//         virtual block<ST, DIMS> set_equal_partition_ith(size_t[DIMS] count, size_t[DIMS] idx) = 0;

//         // fixed partition
//         virtual std::vector<block<ST, DIMS>> set_fixed_partition(ST[DIMS] size) = 0;
//         virtual block<ST, DIMS> set_equal_partition_ith(ST[DIMS] size, size_t[DIMS] idx) = 0;
        
// };



// template <typename ST>
// using block1D = block<ST, 1>;

// template <typename ST>
// class block<ST, 2> {
//     public:
//         ST[2] offset;
//         ST[2] size;

//         virtual bool inside(ST[2] pt) {
//             bool out = true;
//             for (size_t i = 0; i < 2; ++i) {
//                 out &= (pt[i] >= offset[i]) && (pt[i] < (offset[i] + size[i])); 
//             }
//         };

//         // overlap?
//         virtual bool overlap(block<ST, 2> const & rhs) = 0;

//         virtual block<ST, 2>& set_union(block<ST, 2> const & rhs) = 0;
//         // OR operator. union
//         inline block<ST, 2>& operator|=(block<ST, 2> const & rhs) { 
//             return this->set_union(rhs); 
//         };
//         // passing lhs by value helps optimize chained a+b+c.  then return results by value (move contructor, elision)
//         friend block<ST, 2> operator|(block<ST, 2> lhs, block<ST, 2> const & rhs) {
//             return lhs.set_union(lhs);
//         }; 

//         virtual block<ST, 2>& set_intersection(block<ST, 2> const & rhs) = 0;
//         // AND operator.  intersection
//         inline block<ST, 2>& operator&=(block<ST, 2> const & rhs) {
//             return this->set_intersection(rhs);
//         };
//         // passing lhs by value helps optimize chained a+b+c
//         friend block<ST, 2> operator&(block<ST, 2> lhs, block<ST, 2> const & rhs) {
//             return lhs->set_intersection(rhs);
//         }; 

//         virtual block<ST, 2>& set_difference(block<ST, 2> const & rhs) = 0;
//         // set difference operator. 
//         inline block<ST, 2>& operator-=(block<ST, 2> const & rhs) {
//             return this->set_difference(rhs);
//         };
//         // passing lhs by value helps optimize chained a+b+c
//         friend block<ST, 2> operator-(block<ST, 2> lhs, block<ST, 2> const & rhs) {
//             return lhs.set_difference(rhs);
//         }; 

//         // equal partition
//         virtual std::vector<block<ST, 2>> set_equal_partition(size_t[2] count) = 0;
//         virtual block<ST, 2> set_equal_partition_ith(size_t[2] count, size_t[2] idx) = 0;

//         // fixed partition
//         virtual std::vector<block<ST, 2>> set_fixed_partition(ST[2] size) = 0;
//         virtual block<ST, 2> set_equal_partition_ith(ST[2] size, size_t[2] idx) = 0;
        
// };


// template <typename ST>
// using block2D = block<ST, 2>;

template <typename ST>
struct partition {
    static_assert( ::std::is_integral<ST>::value, "Partition: supports only integral template parameter" );

    ST offset;
    ST size;

    partition() = default;
    partition(ST const & _off, ST const & _size) :
        offset(_off), size(_size) {}
    partition(partition const & other) = default;
    partition(partition && other) = default;
    partition& operator=(partition const & other) = default;
    partition& operator=(partition && other) = default;

    void print(const char * prefix) {
        FMT_PRINT_RT("{} Partition: offset: {}, size: {}\n", prefix, offset, size);
    }  
#ifdef USE_MPI
    static partition make_partition(ST const & size, MPI_Comm comm = MPI_COMM_WORLD) {
        int procs;
        int rank;
        MPI_Comm_size(comm, &procs);
        MPI_Comm_rank(comm, &rank);

        ST offset = 0;
        utils::mpi::datatype<ST> dt;
        MPI_Exscan(&size, &offset, 1, dt.value, MPI_SUM, comm);
        return partition(offset, size);
    }
#else
    static partition make_partition(ST const & size) {
        return partition(0, size);
    }
#endif

};

#ifdef USE_MPI
namespace mpi {

template <typename ST>
struct datatype<partition<ST>, false> {
    datatype() {
        utils::mpi::datatype<ST> dt1;
        MPI_Datatype type[2] = {
            dt1.value,
            dt1.value
        };
        int blocklen[2] = {1, 1};
        partition<ST> test {};
        MPI_Aint disp[2]; 
        disp[0] = reinterpret_cast<unsigned char *>(&test.offset) - reinterpret_cast<unsigned char *>(&test);
        disp[1] = reinterpret_cast<unsigned char *>(&test.size) - reinterpret_cast<unsigned char *>(&test);
        MPI_Type_create_struct(2, blocklen, disp, type, &value);
        MPI_Type_commit(&value);
    }
    ~datatype() {
        MPI_Type_free(&value);
    }
    MPI_Datatype value;
};

}
#endif


template <int Strategy>
class partitioner1D;

template <>
class partitioner1D<PARTITION_EQUAL> {

    public:
        // ------------ block partition a 1D range.
        template <typename ST>
        inline std::vector<partition<ST>> divide(ST const & total, ST const & parts) const {
            return divide(partition<ST>(0, total), parts);
        }

        // partitions the given partition.  advance the offset, but does minimal with the parts.
        template <typename ST>
        inline std::vector<partition<ST>> divide(partition<ST> const & src, ST const & parts) const {
            std::vector<partition<ST>> partitions;

            if (parts < 1) return partitions;
            partitions.reserve(parts);
            if (parts == 1) {
                partitions.emplace_back(src);                
                return partitions;
            }

            // more than 1 part.
            ST block = src.size / parts;
            ST remainder = src.size - block * parts;

            ST offset = src.offset;
            ST i = 0;
            block += (remainder > 0);
            for (; i < remainder; ++i) {
                partitions.emplace_back(offset, block);
                offset += block;
            }
            block -= (remainder > 0);
            for (; i < parts; ++i) {
                partitions.emplace_back(offset, block);
                offset += block;
            }

            return partitions;
        }

        // block partition a 1D range and return the target partition
        template <typename ST>
        inline partition<ST> get_partition(ST const & total, int const & parts, int const & idx) const {
            return get_partition(partition<ST>(0, total), parts, idx);
        }

        template <typename ST>
        inline partition<ST> get_partition(partition<ST> const & src, int const & parts, int const & idx) const {
            if ((idx >= parts) || (idx < 0) || (parts < 1)) return partition<ST>(0, 0);
            if (parts == 1) return src;

            ST block = src.size / parts;
            ST remainder = src.size - block * parts;

            return partition<ST>(src.offset + block * idx + (static_cast<ST>(idx) < remainder ? idx : remainder), 
                block + (static_cast<ST>(idx) < remainder ? 1 : 0));
        }

};

template <>
class partitioner1D<PARTITION_FIXED> {

    public:

        // ---------------- fixed size partition a 1D range.
        template <typename ST>
        inline std::vector<partition<ST>> divide(ST const & total, ST const & block) const {
            return divide(partition<ST>(0, total, 0), block);
        }
        template <typename ST>
        inline std::vector<partition<ST>> divide(partition<ST> const & src, ST const & block) const {
            std::vector<partition<ST>> partitions;
            if (block <= 0) return partitions;
            if (block >= src.size) {
                partitions.emplace_back(src);
                return partitions;
            }

            typename partition<ST>::id_type parts = (src.size + block - 1) / block;
            partitions.reserve(parts);

            typename partition<ST>::id_type i = 0;
            ST offset = src.offset;
            for (; i < parts; ++i) {
                partitions.emplace_back(offset, block); 
                offset += block;
            }
            // limit size of the last block.
            partitions[parts - 1].size = src.offset + src.size - partitions[parts - 1].offset;

            return partitions;
        }

        // block partition a 1D range and return the target partition
        template <typename ST>
        inline partition<ST> get_partition(ST const & total, ST const & block, int const & idx) const {
            return get_partition(partition<ST>(0, total, 0), block, idx);
        }
        template <typename ST>
        inline partition<ST> get_partition(partition<ST> const & src, ST const & block, int const & idx) const {
            if ((block <= 0) || (idx < 0)) return partition<ST>(0, 0); 
            typename partition<ST>::id_type parts = (src.size + block - 1) / block;  // number of parts
            if (idx >= parts) return partition<ST>(0, 0); 
            if (parts == 1) return src;

            ST offset = block * idx;   // offset from src.offset

            return partition<ST>( offset + src.offset, 
                std::min(block, src.size - offset));
        }
};

template <typename ST>
struct partition2D {
    static_assert( ::std::is_integral<ST>::value, "Partition: supports only integral template parameter" );

    using part1D_type = partition<ST>;

    part1D_type r;
    part1D_type c;

    partition2D() = default;
    partition2D(part1D_type const & _r, part1D_type const & _c) :
        r(_r), c(_c) {}
    partition2D(partition2D const & other) = default;
    partition2D(partition2D && other) = default;
    partition2D& operator=(partition2D const & other) = default;
    partition2D& operator=(partition2D && other) = default;


    void print(const char * prefix) {
        char pre[1024];
        strcpy(pre, prefix);
        strcpy(pre + strlen(prefix), " ROW ");
        r.print(pre);
        strcpy(pre + strlen(prefix), " COL ");
        c.print(pre);
        FMT_PRINT_RT("{} Partition2D\n", prefix);
    }
};

#ifdef USE_MPI
namespace mpi {

template <typename ST>
struct datatype<partition2D<ST>, false> {
    datatype() {
        utils::mpi::datatype<typename partition2D<ST>::part1D_type> dt1;
        MPI_Datatype type[2] = {
            dt1.value,
            dt1.value
        };
        int blocklen[2] = {1, 1};
        partition2D<ST> test {};
        MPI_Aint disp[2];
        disp[0] = reinterpret_cast<unsigned char *>(&test.r) - reinterpret_cast<unsigned char *>(&test);
        disp[1] = reinterpret_cast<unsigned char *>(&test.c) - reinterpret_cast<unsigned char *>(&test);
        MPI_Type_create_struct(2, blocklen, disp, type, &value);
        MPI_Type_commit(&value);
    }
    ~datatype() {
        MPI_Type_free(&value);
    }
    MPI_Datatype value;
};

}
#endif



template <int Strategy>
class partitioner2D;

template <>
class partitioner2D<PARTITION_EQUAL> {

    protected:
        partitioner1D<PARTITION_EQUAL> r_partitioner;
        partitioner1D<PARTITION_EQUAL> c_partitioner;

    public:
        // -------------- block partition a 2D range, rectilinear.
        template <typename ST>
        inline std::vector<partition2D<ST>> divide(ST const & rows, ST const & cols,
            int const & row_parts, int const & col_parts) const {
            return divide(partition2D<ST>(
                        partition<ST>(0, rows, 0),
                        partition<ST>(0, cols, 0)), 
                        row_parts, 
                        col_parts);
        }
        template <typename ST>
        inline std::vector<partition2D<ST>> divide_rows(ST const & rows, ST const & cols,
            int const & row_parts) const {
            return divide(rows, cols,row_parts, 1);
        }
        template <typename ST>
        inline std::vector<partition2D<ST>> divide_columns(ST const & rows, ST const & cols,
            int const & col_parts) const {
            return divide(rows, cols, 1, col_parts);
        }

        template <typename ST>
        inline std::vector<partition2D<ST>> divide(partition2D<ST> const & src,
            int const & row_parts, int const & col_parts) const {

            int parts = row_parts * col_parts;
            std::vector<partition2D<ST>> partitions;

            if (parts < 1) return partitions;
            partitions.reserve(parts);
            if (parts == 1) {
                partitions.emplace_back(src);
                return partitions;
            }

            // first partition horizontal and vertical
            std::vector<partition<ST>> r_partitions = r_partitioner.divide(src.r, row_parts);
            std::vector<partition<ST>> c_partitions = c_partitioner.divide(src.c, col_parts);

            // iterator and make a combined.
            for (int i = 0; i < row_parts; ++i) {
                for (int j = 0; j < col_parts; ++j) {
                    partitions.emplace_back(
                        r_partitions[i],
                        c_partitions[j]);
                }
            }
            return partitions;
        }

        template <typename ST>
        inline std::vector<partition2D<ST>> divide_rows(partition2D<ST> const & src,
            int const & row_parts) const {
                return divide(src, row_parts, 1);
        }

        template <typename ST>
        inline std::vector<partition2D<ST>> divide_columns(partition2D<ST> const & src,
            int const & col_parts) const {
                return divide(src, 1, col_parts);
        }

        // ----------------- block partition a 2D range rectilinear, and return the target partition
        template <typename ST>
        inline partition2D<ST> get_partition(ST const & rows, ST const & cols,
            int const & row_parts, int const & col_parts,
            int const & row_idx, int const & col_idx) const {
            return get_partition(partition2D<ST>(
                        partition<ST>(0, rows, 0),
                        partition<ST>(0, cols, 0)), 
                        row_parts, col_parts, row_idx, col_idx);
            }

        template <typename ST>
        inline partition2D<ST> get_partition(partition2D<ST> const & src,
            int const & row_parts, int const & col_parts,
            int const & row_idx, int const & col_idx) const {
            if ((row_idx >= row_parts) || (row_idx < 0) || (col_idx >= col_parts) || (col_idx < 0)) return partition<ST>(0, 0);
            int parts = row_parts * col_parts;
            if (parts < 1) return partition<ST>(0, 0);
            if (parts == 1) return src;

            // first partition horizontal and vertical
            partition<ST> r_partition = r_partitioner.get_partition(src.r, row_parts, row_idx);
            partition<ST> c_partition = c_partitioner.get_partition(src.c, col_parts, col_idx);

            return partition2D<ST>(r_partition,
                                    c_partition);
        }

};

template <>
class partitioner2D<PARTITION_FIXED> {
    protected:
        partitioner1D<PARTITION_FIXED> r_partitioner;
        partitioner1D<PARTITION_FIXED> c_partitioner;

    public:

        // ---------------- fixed size partition a 2D range, rectilinear..
        template <typename ST>
        inline std::vector<partition2D<ST>> divide(ST const & rows, ST const & cols,
            ST const & row_block, ST const & col_block) const {
            return divide(partition2D<ST>(
                        partition<ST>(0, rows, 0),
                        partition<ST>(0, cols, 0)), row_block, col_block);
        }
        template <typename ST>
        inline std::vector<partition2D<ST>> divide_rows(ST const & rows, ST const & cols,
            ST const & row_block) const {
            return divide(partition2D<ST>(
                        partition<ST>(0, rows, 0),
                        partition<ST>(0, cols, 0)), row_block, cols);
        }
        template <typename ST>
        inline std::vector<partition2D<ST>> divide_columns(ST const & rows, ST const & cols,
            ST const & col_block) const {
            return divide(partition2D<ST>(
                        partition<ST>(0, rows, 0),
                        partition<ST>(0, cols, 0)), rows, col_block);
        }
        template <typename ST>
        inline std::vector<partition2D<ST>> divide(partition2D<ST> const & src,
            ST const & row_block, ST const & col_block) const {
            std::vector<partition2D<ST>> partitions;
            if ((row_block <= 0) || (col_block <= 0)) return partitions;
            if ((row_block >= src.r.size ) && (col_block >= src.c.size)) {
                partitions.emplace_back(src);
                return partitions;
            }

            // first partition horizontal and vertical
            std::vector<partition<ST>> r_partitions = r_partitioner.divide(src.r, row_block);
            std::vector<partition<ST>> c_partitions = c_partitioner.divide(src.c, col_block);

            size_t row_parts = r_partitions.size();
            size_t col_parts = c_partitions.size();
            partitions.reserve(row_parts * col_parts);

            // iterator and make a combined.
            for (size_t i = 0; i < row_parts; ++i) {
                for (size_t j = 0; j < col_parts; ++j) {
                    partitions.emplace_back(
                        r_partitions[i],
                        c_partitions[j]);
                }
            }

            return partitions;
        }

        template <typename ST>
        inline std::vector<partition2D<ST>> divide_row(partition2D<ST> const & src,
            ST const & row_block) const {
            return divide(src, row_block, src.c.size);
        }
        template <typename ST>
        inline std::vector<partition2D<ST>> divide_col(partition2D<ST> const & src,
            ST const & col_block) const {
            return divide(src, src.r.size, col_block);
        }

        // block partition a 2D range rectilinear, and return the target partition
        template <typename ST>
        inline partition2D<ST> get_partition(ST const & rows, ST const & cols,
            ST const & row_block, ST const & col_block,
            int const & row_idx, int const & col_idx) const {

            return get_partition(partition2D<ST>(
                        partition<ST>(0, rows),
                        partition<ST>(0, cols)), 
                        row_block, col_block, row_idx, col_idx);
        }
                    

        template <typename ST>
        inline partition2D<ST> get_partition(partition2D<ST> const & src,
            ST const & row_block, ST const & col_block,
            int const & row_idx, int const & col_idx) const {

            if ((row_block <= 0) || (row_idx < 0) || (col_block <= 0) || (col_idx < 0)) 
                return partition2D<ST>();

            if ((row_block >= src.r.size) && (col_block >= src.c.size)) return src;

            auto r_parts = (src.r.size + row_block - 1) / row_block;  // number of parts
            auto c_parts = (src.c.size + col_block - 1) / col_block;  // number of parts
            if ((row_idx >= r_parts) || (col_idx >= c_parts)) return partition2D<ST>(); 
            if ((r_parts == 1) && (c_parts == 1)) return src;

            partition<ST> r_partition = r_partitioner.get_partition(src.r, row_block, row_idx);
            partition<ST> c_partition = c_partitioner.get_partition(src.c, col_block, col_idx);

            return partition2D<ST>(
                        r_partition,
                        c_partition);
        }

};


// ------------- additional partition filters.  

// upper triangle filter.  look at offset values to determine if a block is in scope.
class upper_triangle_filter {
    protected:
        int64_t off_diag;   // distance ABOVE diagonal.   0 == diag.  negative dips below diagonal.

    public:
        upper_triangle_filter(int64_t const & dist_from_diag = 0) : off_diag(dist_from_diag) {}

        template <typename ST>
        inline std::vector<partition2D<ST>> filter(std::vector<partition2D<ST>> const & parts,
            partition2D<ST> const & full) const {
            if (parts.size() == 0) return parts;
            
            std::vector<partition2D<ST>> selected;
            selected.reserve(parts.size());

            // keep the original r and c coord, as well as orig cols.  change id to be sequential.
            typename partition2D<ST>::id_type id = 0;
            for (auto part : parts) {
                if (part.c.offset >= (part.r.offset + off_diag)) {
                    selected.push_back(part);
                }
            }
            return selected;
        }

        // if one that will be filtered out, then offset and size = 0.
        template <typename ST>
        partition2D<ST> filter(partition2D<ST> const & part, partition2D<ST> const & full) const {
            if (part.c.offset >= (part.r.offset + off_diag)) return part;
            else return partition2D<ST>();
        }
};


// IMPORTANT ASSUMPTION: WIDTH >= HEIGHT
class banded_diagonal_filter {
    // columns is even or odd, need columns/2 +1 per row.
    // with even number of columns, a little less than full band.
    protected:
        size_t band_width;
        int64_t off_diag;

    public:
        banded_diagonal_filter(
            size_t const & _band_width = std::numeric_limits<size_t>::max(), 
            int64_t const & _dist_from_diag = 0) : 
            band_width(_band_width), off_diag(_dist_from_diag) {}

        template <typename ST>
        inline std::vector<partition2D<ST>> filter(std::vector<partition2D<ST>> const & parts,
            partition2D<ST> const & full) const {
            if (parts.size() == 0) return parts;
            
            std::vector<partition2D<ST>> selected;
            selected.reserve(parts.size());

            // first get min and max column values.
            size_t bw = std::min(this->band_width, static_cast<size_t>(full.c.size / 2 + 1));

            // keep the original r and c coord, as well as orig cols.  change id to be sequential.
            FMT_ROOT_PRINT("Banded Diagnonal Filter for tiles: width = {}, bandwidth = {}\n", full.c.size, bw);
            
            int64_t offdiag = this->off_diag;
            while (offdiag < 0) offdiag += full.c.size;
            offdiag %= full.c.size;   // make sure it's within [0, full.c.size)

            int64_t shifted_c;
            for (auto part : parts) {
                // part.r.offset could be as large as c., same with offdiag, hence add full.c.size 2x.
                shifted_c = static_cast<int64_t>(part.c.offset + full.c.size + full.c.size) - static_cast<int64_t>(part.r.offset + offdiag);
                if (shifted_c >= 0) shifted_c = shifted_c % full.c.size;  // only positive values can handle this.
                // tall matrix:  shifted_c would have negative number.

                if ((shifted_c >= 0) && (shifted_c < bw)) {
                    selected.push_back(part);
                }
            }
            return selected;
        }

        // if one that will be filtered out, then id is -1.
        template <typename ST>
        partition2D<ST> filter(partition2D<ST> const & part, partition2D<ST> const & full) const {

            // first get min and max column values.
            size_t bw = std::min(this->band_width, static_cast<size_t>(full.c.size / 2 + 1));

            // keep the original r and c coord, as well as orig cols.  change id to be sequential.
            FMT_ROOT_PRINT("Banded Diagnonal Filter for tiles: width = {}, bandwidth = {}\n", full.c.size, bw);
            
            int64_t offdiag = this->off_diag;
            while (offdiag < 0) offdiag += full.c.size;
            offdiag %= full.c.size;   // make sure it's within [0, full.c.size)

            int64_t shifted_c = static_cast<int64_t>(part.c.offset + full.c.size + full.c.size) - static_cast<int64_t>(part.r.offset + offdiag);
            if (shifted_c >= 0) shifted_c %= full.c.size;  // only positive values can handle this.
            // tall matrix:  shifted_c would have negative number.

            if ((shifted_c >= 0) && (shifted_c < bw)) return part;
            else return partition2D<ST>();
        }
};




}