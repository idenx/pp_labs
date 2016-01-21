#include "main.hpp"

#include "matrix.hpp"
#include "vector.hpp"
#include "eigenvalue_searcher.hpp"
#include "timer.hpp"

#include <iostream>
#include <fstream>
#include <cstdlib>

bool use_openmp = false;
bool use_mpi = false;
bool use_mpi_timer = false;
bool use_mpi_tracing = false;

int mpi_nodes_n = 0;
int mpi_rank = MPI_RANK_MASTER;

Matrix<matrix_elem_t> g_submatrix;


template <class T>
static Matrix<T> gen_random_matrix(size_t size)
{
//    const T max_value = static_cast<T>(-1);
//rand() % max_value;
    #define GEN_ELEM() m(i, j) = 7 * i + j;

    Matrix<T> m(size, size, 0);

    if (use_openmp) {
        #pragma omp parallel for
        for (int i = 0; i < static_cast<int>(m.rows_n()); ++i)
            for (int j = 0; j < static_cast<int>(m.cols_n()); ++j)
                GEN_ELEM();
    } else {
        for (uint32_t i = 0; i != m.rows_n(); ++i)
            for (uint32_t j = 0; j != m.cols_n(); ++j)
                GEN_ELEM();
    }

    return m;
}

static void init_mpi(int argc, char *argv[])
{
    Timer t("init mpi");
    std::cout << "Use MPI" << std::endl;

    MPI_Init(&argc, &argv);

    if (use_mpi_tracing) {
        char *tracefile = getenv("TVTRACE");
        if (tracefile != NULL) {
            std::cout << "tv tracefile = " << tracefile << std::endl;
            MPI_Pcontrol(TRACEFILES, NULL, tracefile, 0);
        } else {
            MPI_Pcontrol(TRACEFILES, NULL, "trace", 0);
        }

        MPI_Pcontrol(TRACELEVEL, 1, 1, 1);
        MPI_Pcontrol(TRACENODE, 1000000, 1, 1);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &mpi_nodes_n);
    ASSERT(mpi_nodes_n > 1);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    use_mpi_timer = true;
}

// this function is run only by workers
static int do_mpi_work()
{
    ASSERT(mpi_rank != MPI_RANK_MASTER);

    int matrix_size = 0;
    int err;
    {
        Timer t("worker get matrix size");
        err = MPI_Bcast(&matrix_size, 1, MPI_INT, MPI_RANK_MASTER, MPI_COMM_WORLD); // recv matrix size from master
    }
    ASSERT(err == MPI_SUCCESS);
    ASSERT(matrix_size);

    ASSERT(mpi_nodes_n);
    ASSERT(matrix_size >= mpi_nodes_n);
    ASSERT(matrix_size % mpi_nodes_n == 0);

    int rows_per_worker = matrix_size / mpi_nodes_n;
    ASSERT(rows_per_worker);
    int nums_per_worker = rows_per_worker * matrix_size;
    Matrix<matrix_elem_t> submatrix(rows_per_worker, matrix_size, 0);

    {
        Timer t3(string_format("scatter of %zu kb", nums_per_worker * sizeof(matrix_elem_t) / 1024));
        // recv appropriate rows range by each worker
        MpiTracer mt(0);
        err = MPI_Scatter(NULL, 0, MPI_MATRIX_DATA_TYPE, submatrix.data(), nums_per_worker, MPI_MATRIX_DATA_TYPE, MPI_RANK_MASTER, MPI_COMM_WORLD);
        ASSERT(err == MPI_SUCCESS);
    }

    for (;;) {
        int stop;
        MpiTracer mt(1);
        err = MPI_Bcast(&stop, 1, MPI_INT, MPI_RANK_MASTER, MPI_COMM_WORLD);
        ASSERT(err == MPI_SUCCESS);
        if (stop) break;

        std::vector<vector_elem_t> vec(matrix_size);
        err = MPI_Bcast(vec.data(), vec.size(), MPI_VECTOR_DATA_TYPE, MPI_RANK_MASTER, MPI_COMM_WORLD); // send rhs vector to all workers
        ASSERT(err == MPI_SUCCESS);

        ASSERT(use_mpi == false);

        std::vector<vector_elem_t> mul_res;
        // main worker calculation here: it's a multiplication of submatrix by vector
        mul_res = submatrix * vec;
        ASSERT(mul_res.size() == (size_t)rows_per_worker);

        err = MPI_Gather(mul_res.data(), mul_res.size(), MPI_VECTOR_DATA_TYPE, NULL, 0, MPI_VECTOR_DATA_TYPE, MPI_RANK_MASTER, MPI_COMM_WORLD);
        ASSERT(err == MPI_SUCCESS);
    }

    MPI_Finalize();

    return 0;
}

static int do_calculation(int argc, char *argv[], const std::string &mode)
{
    Timer total("total calc");

    if (mode.find("openmp") != std::string::npos) {
        use_openmp = true;
        std::cout << "Use OpenMP with " << omp_get_num_threads() << " threads" << std::endl;
    }

    if (mode.find("mpi") != std::string::npos) {
        use_mpi_tracing = mode.find("trace") != std::string::npos;
        init_mpi(argc, argv);
        if (mpi_rank == MPI_RANK_MASTER) // master
            use_mpi = true;
        else // worker
            return do_mpi_work();
    }

    Matrix<matrix_elem_t> m;
    if (mode.find("input") != std::string::npos) {
        Timer t("read from file");
        std::ifstream in(argv[2]);
        Matrix<matrix_elem_t> m_(in);
        m = m_;
    } else {
        Timer t("gen");
        MpiTracer mt(4);
        size_t size = atoi(argv[2]);
        m = gen_random_matrix<matrix_elem_t>(size);
    }

    EigenvalueSearcher<matrix_elem_t> es;
    std::vector<vector_elem_t> evec;
    vector_elem_t eval;

    {
        Timer t("calc");
        if (use_mpi) {
            // broadcast matrix size
            int matrix_size = m.rows_n();
            {
                Timer t_bcast("bcast matrix size");
                int err = MPI_Bcast(&matrix_size, 1, MPI_INT, MPI_RANK_MASTER, MPI_COMM_WORLD);
                ASSERT(err == MPI_SUCCESS);
            }

            // send appropriate rows range to each worker
            {
                Timer t_scatter_marix("scatter matrix");
                int rows_per_worker = matrix_size / mpi_nodes_n;
                ASSERT(rows_per_worker);
                int nums_per_worker = rows_per_worker * matrix_size;
                g_submatrix = Matrix<matrix_elem_t>(rows_per_worker, m.cols_n(), 0);
                int err = MPI_Scatter(m.data(), nums_per_worker, MPI_MATRIX_DATA_TYPE, g_submatrix.data(),
                                  nums_per_worker, MPI_MATRIX_DATA_TYPE, MPI_RANK_MASTER, MPI_COMM_WORLD);
                ASSERT(err == MPI_SUCCESS);
            }
        }

        eval = es.get_eigenvalue(m, evec);
    }

    if (mode.find("output") != std::string::npos) {
        Timer t("print to file");
        std::ofstream out(argv[3]);
        out << eval << " " << std::endl << evec;
    }

    std::cout << "found eigenvector and max eigenvalue " << eval << " for matrix of size " << m.rows_n() << std::endl;
    if (use_mpi) {
        int stop = 1;
        int err = MPI_Bcast(&stop, 1, MPI_INT, MPI_RANK_MASTER, MPI_COMM_WORLD);
        ASSERT(err == MPI_SUCCESS);

        MPI_Finalize();
    }

    return 0;
}

static int do_generation(int argc, char *argv[])
{
    Timer t("generation");
    std::string fname(argv[2]);
    std::ofstream out(fname.c_str());
    size_t size = atoi(argv[3]);
    Matrix<matrix_elem_t> m = gen_random_matrix<matrix_elem_t>(size);
    std::cout << "generated matrix of size " << size << ", saving it into " << fname << " ..." << std::endl;
    out << m;
    return 0;
}

int main(int argc, char *argv[])
{
    srand(0);

    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " mode ..." << std::endl;
        return EXIT_FAILURE;
    }

    std::string mode(argv[1]);

    if (mode.find("calc") != std::string::npos)
        return do_calculation(argc, argv, mode);

    if (mode == "gen")
        return do_generation(argc, argv);

    return EXIT_FAILURE;
}
