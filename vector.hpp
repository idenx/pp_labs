#ifndef _VECTOR_HPP_
#define _VECTOR_HPP_
#include "main.hpp"
#include <cmath>

template <typename T, typename R>
std::vector<T> operator*(const std::vector<T> &vec, R scalar)
{
    std::vector<T> res_vec(vec);
    for (typename std::vector<T>::iterator i = res_vec.begin(); i != res_vec.end(); ++i)
        *i *= scalar;

    return res_vec;
}

template <typename T, typename R>
std::vector<T> operator/=(std::vector<T> &vec, R scalar)
{
    ASSERT(scalar);
    for (typename std::vector<T>::iterator i = vec.begin(); i != vec.end(); ++i)
        *i /= scalar;

    return vec;
}


template <typename T>
std::vector<T> operator-(const std::vector<T> &a, const std::vector<T> &b)
{
    std::vector<T> res_vec;
    ASSERT(a.size() == b.size());
    res_vec.resize(a.size());

    for (typename std::vector<T>::size_type i = 0; i != a.size(); ++i)
        res_vec[i] = a[i] - b[i];

    return res_vec;
}

template <typename T>
std::ostream& operator<<(std::ostream &os, const std::vector<T> &vec)
{
    for (typename std::vector<T>::const_iterator i = vec.begin(); i != vec.end(); ++i)
        os << *i << " ";

    return os;
}

template <typename T>
T vector_get_max_abs_value(const std::vector<T> &vec)
{
    ASSERT(vec.size());
    T max_value = std::fabs(vec[0]);
    for (typename std::vector<T>::const_iterator i = vec.begin(); i != vec.end(); ++i) {
        if (std::fabs(*i) > max_value)
            max_value = std::fabs(*i);
    }

    return max_value;
}

#endif
