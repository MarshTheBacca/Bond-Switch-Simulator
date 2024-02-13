#include "vecr.h"
#include <iostream>
#include <sstream>
// ##### VECTOR RESERVED SIZE CLASS #####

/**
 * @brief Default constructor for an empty VecR object with a maximum size of 1
 * @typedef T The type of the vector
 */
template <typename T>
VecR<T>::VecR() {
    this->n = 0;
    this->nMax = 1;
    this->v = new T[this->nMax]();
}

/**
 * @brief Construct a new VecR object with a given maximum size
 * @param maxSize The maximum size of the vector
 * @typedef T The type of the vector
 */
template <typename T>
VecR<T>::VecR(int maxSize) {
    if (maxSize < 0)
        throw std::runtime_error("Vector cannot be instantiated with negative size");
    this->n = maxSize;
    this->nMax = maxSize;
    this->v = new T[this->nMax]();
}

/**
 * @brief Construct a new VecR object with a given size and maximum size
 * @param size The size of the vector
 * @param maxSize The maximum size of the vector
 * @typedef T The type of the vector
 */
template <typename T>
VecR<T>::VecR(int size, int maxSize) {
    if (maxSize < 0)
        throw std::runtime_error("Vector cannot be instantiated with negative size");
    this->n = size;
    this->nMax = maxSize;
    this->v = new T[this->nMax]();
}

/**
 * @brief Construct a new VecR object by copying another
 * @param source The source VecR object
 * @typedef T The type of the vector
 */
template <typename T>
VecR<T>::VecR(const VecR &source) {
    this->n = source.n;
    this->nMax = source.nMax;
    this->v = new T[this->nMax]();
    for (int i = 0; i < this->n; ++i)
        this->v[i] = source.v[i];
}

/**
 * @brief Destroy the VecR object
 * @typedef T The type of the vector
 */
template <typename T>
VecR<T>::~VecR() {
    delete[] this->v;
}

// Check for equality of any value type
template <typename T>
inline bool VecR<T>::equals(const T &a, const T &b) const {
    return a == b || abs(a - b) < abs(std::min(a, b)) * std::numeric_limits<T>::epsilon();
}

/**
 * @brief Set the current size of the vector
 * @param size The new size of the vector
 * @typedef T The type of the vector
 */
template <typename T>
void VecR<T>::setSize(int size) {
    if (size > this->nMax) {
        std::stringstream ss;
        ss << "Cannot set vector size larger than reserved size, nMax = " << this->nMax;
        throw std::runtime_error(ss.str());
    }
    this->n = size;
}

/**
 * @brief Reset the maximum size of the vector
 * @param maxSize The new maximum size of the vector, setting additional values to 0
 * @typedef T The type of the vector
 */
template <typename T>
void VecR<T>::resetMaxSize(int maxSize) {
    if (maxSize != this->nMax) {
        auto newV = std::make_unique<T[]>(maxSize);
        std::copy(this->v, this->v + this->n, newV.get());
        this->nMax = maxSize;
        this->v = newV.release();
    }
}

/**
 * @brief Add a value to the end of the vector
 * @param value The value to add
 * @typedef T The type of the vector
 */
template <typename T>
void VecR<T>::addValue(T value) {
    if (this->n == this->nMax) {
        std::stringstream ss;
        ss << "Cannot add to vector to make larger than reserved size, nMax = " << this->nMax;
        throw std::runtime_error(ss.str());
    }
    this->v[this->n] = value;
    ++this->n;
}

/**
 * @brief Delete the first instance of a value from the vector
 * @param value The value to delete
 * @typedef T The type of the vector
 */
template <typename T>
void VecR<T>::delValue(T value) {
    int d;
    for (d = 0; d < this->n; ++d) {
        if (equals(this->v[d], value)) {
            break;
        }
    }
    if (d == this->n) {
        throw std::runtime_error("Cannot delete value as not present in vector");
    }
    std::move(this->v + d + 1, this->v + this->n, this->v + d);
    --this->n;
}

/**
 * @brief Replace a value in the vector with another
 * @param removeValue The value to remove
 * @param addValue The value to add
 * @param replaceAll Whether to replace all instances of the value
 * @typedef T The type of the vector
 * @throws std::runtime_error If the value to replace is not present in the vector
 */
template <typename T>
void VecR<T>::replaceValue(T removeValue, T addValue, bool replaceAll) {
    // Use std::find_if to find the first occurrence of removeValue in the vector
    auto it = std::find_if(this->v, this->v + this->n, [&](const T &value) {
        return equals(value, removeValue);
    });

    // If removeValue is not found in the vector, throw an exception
    if (it == this->v + this->n) {
        std::ostringstream oss;
        oss << "Cannot replace value " << removeValue << " with " << addValue << " as not present in vector: " << *this;
        throw std::runtime_error(oss.str());
    }

    // Replace the found value with addValue
    *it = addValue;

    // If replaceAll is true, continue to find and replace all subsequent occurrences of removeValue
    if (replaceAll) {
        while ((it = std::find_if(it + 1, this->v + this->n, [&](const T &value) {
                    return equals(value, removeValue);
                })) != this->v + this->n) {
            *it = addValue;
        }
    }
}

template <typename T>
bool VecR<T>::tryReplace(int i, int j, int k, T removeValue, T addValue, T betweenValue1, T betweenValue2) {
    if (equals(v[i], removeValue)) {
        if ((equals(v[j], betweenValue1) && equals(v[k], betweenValue2)) ||
            (equals(v[j], betweenValue2) && equals(v[k], betweenValue1))) {
            v[i] = addValue;
            return true;
        }
    }
    return false;
}

/**
 * @brief Replace a value in place of another at a specific point
 * @param removeValue The value to delete
 * @param addValue The value to add
 * @param betweenValue1 The first value between which to swap
 * @param betweenValue2 The second value between which to swap
 * @typedef T The type of the vector
 * @throws std::runtime_error If the sequence is not present in the vector
 */
template <typename T>
void VecR<T>::replaceValue(T removeValue, T addValue, T betweenValue1, T betweenValue2) {
    if (tryReplace(0, n - 1, 1, removeValue, addValue, betweenValue1, betweenValue2) ||
        tryReplace(n - 1, n - 2, 0, removeValue, addValue, betweenValue1, betweenValue2)) {
        return;
    }

    for (int i = 1, j = 0, k = 2; i < n - 1; ++i, ++j, ++k) {
        if (tryReplace(i, j, k, removeValue, addValue, betweenValue1, betweenValue2)) {
            return;
        }
    }
    std::ostringstream oss;
    oss << "Cannot replace value " << removeValue << " with value " << addValue << " between " << betweenValue1 << " and " << betweenValue2 << " as sequence not present in vector: " << *this;
    throw std::runtime_error(oss.str());
}

// Insert value in vector between two others
template <typename T>
void VecR<T>::insertValue(T vInsert, T vBetween0, T vBetween1) {
    if (n == nMax)
        throw std::runtime_error("Cannot insert to vector to make larger than reserved size");
    bool insert = false;
    int insertPos = -1;

    for (int i = 0; i < n; ++i) {
        if ((v[i] == vBetween0 && v[(i + 1) % n] == vBetween1) || (v[i] == vBetween1 && v[(i + 1) % n] == vBetween0)) {
            insertPos = (i + 1) % n;
            insert = true;
        }
        if (insert)
            break;
    }
    if (!insert) {
        std::ostringstream oss;
        oss << "Cannot insert " << vInsert << " between " << vBetween0 << " and " << vBetween1 << " in vector: " << *this;
        throw std::runtime_error(oss.str());
    }
    for (int i = n; i > insertPos; --i)
        v[i] = v[i - 1];
    v[insertPos] = vInsert;
    ++n;
}

// Binary Operators with constant
template <typename T>
void VecR<T>::operator=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] = k;
}

template <typename T>
void VecR<T>::operator+=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] += k;
}

template <typename T>
void VecR<T>::operator-=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] -= k;
}

template <typename T>
void VecR<T>::operator*=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] *= k;
}

template <typename T>
void VecR<T>::operator/=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] /= k;
}

template <typename T>
bool VecR<T>::operator==(const T &k) {
    bool equality = true;
    for (int i = 0; i < this->n; ++i) {
        if (!equals(this->v[i], k))
            equality = false;
    }
    return equality;
}

template <typename T>
bool VecR<T>::operator<(const T &k) {
    bool lt = true;
    for (int i = 0; i < this->n; ++i) {
        if (this->v[i] >= k)
            lt = false;
    }
    return lt;
}

template <typename T>
bool VecR<T>::operator>(const T &k) {
    bool gt = true;
    for (int i = 0; i < this->n; ++i) {
        if (this->v[i] <= k)
            gt = false;
    }
    return gt;
}

template <typename T>
VecR<T> VecR<T>::operator+(const T &k) {
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] + k;
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator-(const T &k) {
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] - k;
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator*(const T &k) {
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] * k;
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator/(const T &k) {
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] / k;
    return vec;
}

// Binary operators with VecR
template <typename T>
void VecR<T>::operator+=(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] += source.v[i];
}

template <typename T>
void VecR<T>::operator-=(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] -= source.v[i];
}

template <typename T>
void VecR<T>::operator*=(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] *= source.v[i];
}

template <typename T>
void VecR<T>::operator/=(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] /= source.v[i];
}

template <typename T>
bool VecR<T>::operator==(const VecR &source) {
    bool equality = true;
    if (this->n != source.n)
        equality = false;
    else {
        for (int i = 0; i < this->n; ++i) {
            if (!equals(this->v[i], source.v[i]))
                equality = false;
        }
    }
    return equality;
}

template <typename T>
VecR<T> &VecR<T>::operator=(const VecR &source) {
    if (this != &source) {
        delete[] this->v;
        this->n = source.n;
        this->nMax = source.nMax;
        this->v = new T[this->nMax]();
        for (int i = 0; i < this->n; ++i)
            this->v[i] = source.v[i];
    }
    return *this;
}

template <typename T>
VecR<T> VecR<T>::operator+(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] + source.v[i];
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator-(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] - source.v[i];
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator*(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] * source.v[i];
    return vec;
}

template <typename T>
VecR<T> VecR<T>::operator/(const VecR &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] / source.v[i];
    return vec;
}

// Unary operators
template <typename T>
VecR<T> VecR<T>::operator-() {
    VecR<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = -this->v[i];
    return vec;
}

/**
 * @brief Check if value is in vector
 * @param value The value to check
 * @return True if value is in vector
 */
template <typename T>
bool VecR<T>::contains(T value) {
    for (int i = 0; i < this->n; ++i) {
        if (this->v[i] == value)
            return true;
    }
    return false;
}

template <typename T>
std::string VecR<T>::toString() {
    std::ostringstream os;
    os << "VecR<" << typeid(T).name() << "> of size " << this->n << ": [";
    for (const T &value : *this) {
        os << value << ' ';
    }
    os << "]";
    return os.str();
}
// Subscript operator
template <typename T>
T &VecR<T>::operator[](int i) {
    if (n <= i) {
        std::ostringstream oss;
        oss << "Vector subscript out of bounds: " << i << " >= " << n << " for vector: " << *this;
        throw std::runtime_error(oss.str());
    }
    return v[i];
}

// Subscript operator
template <typename T>
T &VecR<T>::operator[](int i) const {
    if (n <= i) {
        std::ostringstream oss;
        oss << "Vector subscript out of bounds: " << i << " >= " << n << " for vector: " << *this;
        throw std::runtime_error(oss.str());
    }
    return v[i];
}

/**
 * @brief Stream operator for VecR
 * @param os Output stream
 * @param vec Vector to output
 * @typedef U Type of vector
 * @return Output stream
 */
template <typename U>
std::ostream &operator<<(std::ostream &os, VecR<U> &vec) {
    os << "VecR<" << typeid(U).name() << "> of size " << vec.n << " and max size " << vec.nMax << ": [";
    for (const U &value : vec) {
        os << value << ' ';
    }
    os << "]";
    return os;
}

/**
 * @brief Const version of stream operator for VecR
 * @param os Output stream
 * @param vec Vector to output
 * @typedef U Type of vector
 * @return Output stream
 */
template <typename U>
std::ostream &operator<<(std::ostream &os, const VecR<U> &vec) {
    os << "VecR<" << typeid(U).name() << "> of size " << vec.n << "and max size " << vec.nMax << ": [";
    for (const U &value : vec) {
        os << value << ' ';
    }
    os << "]";
    return os;
}