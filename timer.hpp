#ifndef _TIMER_HPP_
#define _TIMER_HPP_

#include <sys/time.h>
#include "main.hpp"

class Timer {
private:
    struct timeval start_ts_;
    double mpi_start_ts_;
    std::string name_;

public:
    Timer(const std::string &name) : name_(name) {
        if (use_mpi_timer) {
            mpi_start_ts_ = MPI_Wtime();
        } else {
            mpi_start_ts_ = 0;
            gettimeofday(&start_ts_, NULL);
        }
    }
    ~Timer() {
        uint64_t elapsed_ms;
        if (mpi_start_ts_) {
            elapsed_ms = static_cast<uint64_t>((MPI_Wtime() - mpi_start_ts_) * 1000);
        } else {
            struct timeval ts;
            gettimeofday(&ts, NULL);
            elapsed_ms = (ts.tv_sec - start_ts_.tv_sec) * 1000 + (ts.tv_usec - start_ts_.tv_usec) / 1000;
        }
        std::cout << string_format("[%d] Timer %s: elapsed %u ms\n", mpi_rank, name_.c_str(), (uint32_t)elapsed_ms);
    }
};
#endif
