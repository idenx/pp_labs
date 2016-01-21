#ifndef _EIGENVALUE_SEARCHER_HPP_
#define _EIGENVALUE_SEARCHER_HPP_

#include "main.hpp"
#include "eigenvalue_searcher_test.hpp"
#include "vector.hpp"

#include <cmath>
#include <iostream>

template <typename T>
class EigenvalueSearcher {
public:
    template <typename V>
    V get_eigenvalue(const Matrix<T> &matrix, std::vector<V> &eigenvector) {
        // calculate
        eigenvector.resize(matrix.cols_n(), 1);
        V prev_eigenvalue = 0;
        V cur_eigenvalue = 0;
        int index_no = 0;

        V cur_tolerance = 1e9;
        ASSERT(cur_tolerance > tolerance);

        int iter_no = 0;
        while (cur_tolerance > tolerance) {
            V elem = eigenvector[index_no];
            eigenvector = matrix * eigenvector;
            prev_eigenvalue = cur_eigenvalue;
            cur_eigenvalue = eigenvector[index_no] / elem;

            V max_value_of_eigenvector = vector_get_max_abs_value(eigenvector);
            if (max_value_of_eigenvector == 0)
                std::cerr << eigenvector << std::endl;
            ASSERT(max_value_of_eigenvector);
            eigenvector /= max_value_of_eigenvector;

            cur_tolerance = std::fabs(cur_eigenvalue - prev_eigenvalue);
            ++iter_no;
            if (iter_no % 1000 == 0)
                std::cout << "#" << iter_no + 1 << " cur_tolerance: " << cur_tolerance << ", eigenvalue: " << cur_eigenvalue << ", eigenvector after mul: " << eigenvector << std::endl;
        }
        std::cout << "found eigenvalue for " << iter_no << " iterations, cur tolerance is " << cur_tolerance << std::endl;

        // test calculations
        test_eigenvalue_searcher(matrix, cur_eigenvalue, eigenvector);

        return cur_eigenvalue;
    }
};

#endif
