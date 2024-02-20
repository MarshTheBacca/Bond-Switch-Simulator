// Vectors class with fixed size
#ifndef NL_VECF_H
#define NL_VECF_H

#include <cmath>
#include <iostream>
#include <limits>
#include <spdlog/spdlog.h>
#include <sstream>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

// Vector class
template <typename T>
class VecF {

  public:
    // Data members
    int n; // number of values
    T *v;  // values

    // Constructors, copy constructor, destructor
    VecF();
    explicit VecF(int size);
    VecF(const VecF &source);
    ~VecF();

    T *begin() { return &v[0]; }
    T *end() { return &v[n]; }
    const T *begin() const { return &v[0]; }
    const T *end() const { return &v[n]; }

    // Member functions
    bool equals(const T &a, const T &b) const; // check for equality

    // Subscript operator
    T &operator[](int i);
    T &operator[](int i) const;
    // Binary operators with constant
    void operator=(const T &k);
    void operator+=(const T &k);
    void operator-=(const T &k);
    void operator*=(const T &k);
    void operator/=(const T &k);
    bool operator==(const T &k);
    bool operator<(const T &k);
    bool operator>(const T &k);
    VecF operator+(const T &k);
    VecF operator-(const T &k);
    VecF operator*(const T &k);
    VecF operator/(const T &k);
    // Binary operators with vector
    void operator+=(const VecF &source);
    void operator-=(const VecF &source);
    void operator*=(const VecF &source);
    void operator/=(const VecF &source);
    bool operator==(const VecF &source);
    VecF &operator=(const VecF &source);
    VecF operator+(const VecF &source);
    VecF operator-(const VecF &source);
    VecF operator*(const VecF &source);
    VecF operator/(const VecF &source);
    // Unary operators
    VecF operator-();

    bool contains(T value);
    std::string toString();

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, VecF<U> &vec);

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const VecF<U> &vec);

    /**
     * @brief Converts VecF to std::vector
     * @return std::vector of type T
     */
    operator std::vector<T>() const {
        std::vector<T> vec(this->n);
        for (int i = 0; i < this->n; ++i)
            vec[i] = this->v[i];
        return vec;
    }
};

#include "vecf.tpp"
#endif // NL_VECF_H
