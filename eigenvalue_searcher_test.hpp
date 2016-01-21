#ifndef _EIGENVALUE_SEARCHER_TEST_HPP_
#define _EIGENVALUE_SEARCHER_TEST_HPP_

#include "main.hpp"
#include "matrix.hpp"
#include "eigenvalue_searcher.hpp"
#include "vector.hpp"
#include <cmath>

template <typename V>
static void check_vector_components_le(const std::vector<V> &v, V cmp_val)
{
    for (typename std::vector<V>::const_iterator i = v.begin(); i != v.end(); ++i) {
        if (std::abs(*i) > cmp_val) {
            std::cout << "vector is bigger than cmp_val " << cmp_val << ":" << v << std::endl;
            abort();
        }
    }
}

template <typename V, typename M>
void test_eigenvalue_searcher(const Matrix<M> &m, V eigenvalue, const std::vector<V> &eigenvector)
{
    ASSERT(eigenvalue != 0);

    std::vector<V> mul_res = m * eigenvector;
    std::vector<V> mul_eigenvector = eigenvector * eigenvalue;

    std::vector<V> diff = mul_res - mul_eigenvector;
    check_vector_components_le(diff, tolerance);
}

#endif
