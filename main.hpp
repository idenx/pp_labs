#ifndef _MAIN_HPP_
#define _MAIN_HPP_

#include <mpi.h>
#include <pcontrol.h>

#include <assert.h>
#include <iostream>

extern bool use_openmp;
extern bool use_mpi;
extern bool use_mpi_timer;
extern bool use_mpi_tracing;

extern int mpi_nodes_n;
extern int mpi_rank;

const int MPI_RANK_MASTER = 0;

typedef double vector_elem_t;
#define MPI_VECTOR_DATA_TYPE MPI_DOUBLE

typedef uint8_t matrix_elem_t;
#define MPI_MATRIX_DATA_TYPE MPI_BYTE

const vector_elem_t tolerance = 1e-5;

#define string_format(fmt_, args_...) ({ \
    char buf_[4096]; \
    snprintf(buf_, sizeof(buf_), fmt_, ##args_); \
    buf_; \
})


#define ASSERT(c_) ({ \
    if (!(c_)) { \
        std::cerr << "Condition " << #c_ << " failed: " << std::endl; \
        exit(EXIT_FAILURE); \
    } \
})

#define ASSERT_EX(c_, fmt_, args_...) ({ \
    if (!(c_)) { \
        std::cerr << "Condition " << #c_ << " failed: " << string_format(fmt_, ##args_) << std::endl; \
        exit(EXIT_FAILURE); \
    } \
})

#define MPI_TRACE(args_...) ({ if (use_mpi_tracing) MPI_Pcontrol(TRACEEVENT, ##args_); })

class MpiTracer {
private:
    int color_;

public:
    MpiTracer(int color) : color_(color) { MPI_TRACE("entry", color_, 0, ""); }
    ~MpiTracer() { MPI_TRACE("exit", color_, 0, ""); }
};

#endif
