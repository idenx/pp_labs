#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#include <mpi.h>
#include <vector>
#include <ostream>

#define ELEM(i_, j_) mat[i_ * cols + j_]

template <typename T> class Matrix;

extern Matrix<matrix_elem_t> g_submatrix;

template <typename T>
class Matrix {
    private:
        T *mat;
        unsigned rows;
        unsigned cols;

        template <typename V>
        void matrix_mul_vector(std::vector<V> &result, const std::vector<V>& rhs) const {
            ASSERT(rhs.size() == cols);
            ASSERT(result.size() == rows);

            if (use_openmp) {
                #pragma omp parallel for
                for (int i = 0; i < static_cast<int>(rows); i++) {
                    for (int j = 0; j < static_cast<int>(cols); j++) {
                        result[i] += static_cast<V>(ELEM(i, j)) * rhs[j];
                    }
                }
                return;
            }

            for (unsigned i=0; i<rows; i++) {
                for (unsigned j=0; j<cols; j++) {
                    result[i] += static_cast<V>(ELEM(i, j)) * rhs[j];
                }
            }
        }

        template<typename V>
        void matrix_mul_vector_mpi(std::vector<V> &result, const std::vector<V>& rhs) const {
            ASSERT_EX(rows == cols, "rows = %u, cols = %u", rows, cols);
            int matrix_size = rows;
            ASSERT(mpi_rank == MPI_RANK_MASTER); // workers execute another code

            ASSERT(sizeof(T) == sizeof(matrix_elem_t)); // XXX: hack, detemplating here
            ASSERT(mpi_nodes_n);
            ASSERT(matrix_size >= mpi_nodes_n);
            ASSERT(matrix_size % mpi_nodes_n == 0);

            int rows_per_worker = matrix_size / mpi_nodes_n;
            ASSERT(rows_per_worker);

            MpiTracer mt(2);
            int stop = 0;
            int err = MPI_Bcast(&stop, 1, MPI_INT, MPI_RANK_MASTER, MPI_COMM_WORLD);
            ASSERT(err == MPI_SUCCESS);

            // broadcast rhs vector to all workers
            MPI_Bcast(static_cast<void *>(const_cast<V *>(rhs.data())), rhs.size(), MPI_VECTOR_DATA_TYPE, MPI_RANK_MASTER, MPI_COMM_WORLD);

            // master's part of calculation (use global g_submatrix here)
            std::vector<V> mul_res(rows_per_worker);
            ASSERT((int)mul_res.size() == rows_per_worker);
            g_submatrix.matrix_mul_vector(mul_res, rhs);

            // gather all submultiplications from workers
            err = MPI_Gather(mul_res.data(), mul_res.size(), MPI_VECTOR_DATA_TYPE, result.data(),
                             rows_per_worker, MPI_VECTOR_DATA_TYPE, MPI_RANK_MASTER, MPI_COMM_WORLD);
            ASSERT(err == MPI_SUCCESS);
        }


    public:
        Matrix() : mat(NULL), rows(0), cols(0) {};
        Matrix(std::istream &is);
        Matrix(unsigned _rows, unsigned _cols, const T& _initial);
        Matrix(const Matrix<T>& rhs);
        virtual ~Matrix();

        // Operator overloading, for "standard" mathematical matrix operations
        Matrix<T>& operator=(const Matrix<T>& rhs);

        template<typename V>
        std::vector<V> operator*(const std::vector<V>& rhs) const {
            ASSERT(rhs.size() == cols);
            std::vector<V> result(rows, 0);

            if (use_mpi) {
                matrix_mul_vector_mpi(result, rhs);
                return result;
            }

            matrix_mul_vector(result, rhs);

            return result;
        }


        // Access the individual elements
        T& operator()(const unsigned& row, const unsigned& col);
        const T& operator()(const unsigned& row, const unsigned& col) const;

        // Access the row and column sizes
        unsigned rows_n() const;
        unsigned cols_n() const;

        T *data() const { return mat; }
        T *data() { return mat; }
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const Matrix<T> &m)
{
    os << m.rows_n() << " " << m.cols_n() << std::endl;
    for (uint32_t i = 0; i != m.rows_n(); ++i) {
        for (uint32_t j = 0; j != m.cols_n(); ++j) {
            os << m(i, j) << " ";
        }
        os << std::endl;
    }

    return os;
}

#include "matrix.cpp"

#endif
