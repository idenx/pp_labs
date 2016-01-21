#include <omp.h>
#include "main.hpp"
#include "timer.hpp"

// Parameter Constructor
template<typename T>
Matrix<T>::Matrix(unsigned _rows, unsigned _cols, const T& _initial) {
    rows = _rows;
    cols = _cols;
    mat = new T[rows * cols];
    for (unsigned i = 0; i < rows; ++i) {
        for (unsigned j = 0; j < cols; ++j) {
            ELEM(i, j) = _initial;
        }
    }
}

template <typename T>
Matrix<T>::Matrix(std::istream &is)
{
    std::string rows_str, cols_str;
    is >> rows_str >> cols_str;
    rows = atoi(rows_str.c_str());
    cols = atoi(cols_str.c_str());

    mat = new T[rows * cols];
    for (unsigned i = 0; i < rows; ++i) {
        for (unsigned j = 0; j < cols; ++j) {
            is >> ELEM(i, j);
        }
    }
}

// Copy Constructor
template<typename T>
Matrix<T>::Matrix(const Matrix<T>& rhs) {
    mat = new T[rows * cols];
    memcpy(mat, rhs.mat, rows * cols * sizeof(T));
    rows = rhs.rows_n();
    cols = rhs.cols_n();
}

// (Virtual) Destructor
template<typename T>
Matrix<T>::~Matrix() {
    if (mat)
        delete[] mat;
}

// Assignment Operator
template<typename T>
Matrix<T>& Matrix<T>::operator=(const Matrix<T>& rhs) {
    if (&rhs == this)
        return *this;

    if (mat)
        delete[] mat;

    rows = rhs.rows_n();
    cols = rhs.cols_n();
    mat = new T[rows * cols];
    memcpy(mat, rhs.mat, rows * cols * sizeof(T));

    return *this;
}

// Access the individual elements
template<typename T>
T& Matrix<T>::operator()(const unsigned& row, const unsigned& col) {
    return ELEM(row, col);
}

// Access the individual elements (const)
template<typename T>
const T& Matrix<T>::operator()(const unsigned& row, const unsigned& col) const {
    return ELEM(row, col);
}

// Get the number of rows of the matrix
template<typename T>
unsigned Matrix<T>::rows_n() const {
    return this->rows;
}

// Get the number of columns of the matrix
template<typename T>
unsigned Matrix<T>::cols_n() const {
    return this->cols;
}
