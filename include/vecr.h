// Vector class with reserved size
#ifndef NL_VECR_H
#define NL_VECR_H

#include <cmath>
#include <iostream>
#include <limits>
#include <spdlog/spdlog.h>
#include <sstream>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

/* Vector class with reserved size
 * Maximum size defines limit on the number of values
 * Current size determines accessible number of values */
template <typename T>
class VecR {

  public:
    // Data members
    int n;    // number of values
    int nMax; // maximum number of values
    T *v;     // values

    // Constructors, copy constructor, destructor
    VecR();
    explicit VecR(int maxSize);
    VecR(int size, int maxSize);
    VecR(const VecR &source);
    ~VecR();

    // Member functions
    bool equals(const T &a, const T &b) const;                                      // check for equality
    void setSize(int size);                                                         // set current size
    void resetMaxSize(int maxSize);                                                 // reset maximum size
    void addValue(T value);                                                         // add value to end
    void delValue(T value);                                                         // remove first instance of value from vector
    void replaceValue(T vDel, T vAdd, bool swapAll = true);                         // swap a value in place of another
    bool tryReplace(int i, int j, int k, T vDel, T vAdd, T vBetween0, T vBetween1); // Helper function for replaceValue
    void replaceValue(T vDel, T vAdd, T vBetween0, T vBetween1);                    // swap value in place of another at specific point
    void insertValue(T vInsert, T vBetween0, T vBetween1);                          // insert a value between two others

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
    VecR operator+(const T &k);
    VecR operator-(const T &k);
    VecR operator*(const T &k);
    VecR operator/(const T &k);
    // Binary operators with vector
    void operator+=(const VecR &source);
    void operator-=(const VecR &source);
    void operator*=(const VecR &source);
    void operator/=(const VecR &source);
    bool operator==(const VecR &source);
    VecR &operator=(const VecR &source);
    VecR operator+(const VecR &source);
    VecR operator-(const VecR &source);
    VecR operator*(const VecR &source);
    VecR operator/(const VecR &source);
    // Unary operators
    VecR operator-();

    T *begin() { return v; }
    const T *begin() const { return v; } // const version

    T *end() { return v + n; }
    const T *end() const { return v + n; } // const version

    T &back() {
        if (n > 0) {
            return v[n - 1];
        } else {
            throw std::out_of_range("VecR is empty.");
        }
    }

    const T &back() const {
        if (n > 0) {
            return v[n - 1];
        } else {
            throw std::out_of_range("VecR is empty.");
        }
    }

    bool contains(T value);
    std::string toString();

    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, VecR<U> &vec);
    template <typename U>
    friend std::ostream &operator<<(std::ostream &os, const VecR<U> &vec);

    /**
     * @brief Converts VecR to std::vector
     * @return std::vector of type T
     */
    operator std::vector<T>() const {
        std::vector<T> vec(this->n);
        for (int i = 0; i < this->n; ++i)
            vec[i] = this->v[i];
        return vec;
    }

    /**
     * @brief Converts VecR to std::vector
     * @return std::vector of type T
     */
    operator std::vector<T>() {
        std::vector<T> vec(this->n);
        for (int i = 0; i < this->n; ++i)
            vec[i] = this->v[i];
        return vec;
    }
};

#include "vecr.tpp"

#endif // NL_VECR_H
