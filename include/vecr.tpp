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

// Swaps value from vector
template <typename T>
void VecR<T>::swapValue(T vDel, T vAdd, bool swapAll) {
    bool swap = false;
    if (swapAll) {
        for (int i = 0; i < this->n; ++i) {
            if (equals(this->v[i], vDel)) {
                this->v[i] = vAdd;
                swap = true;
            }
        }
    } else {
        for (int i = 0; i < this->n; ++i) {
            if (equals(this->v[i], vDel)) {
                this->v[i] = vAdd;
                swap = true;
                break;
            }
        }
    }
    if (!swap)
        throw std::runtime_error("Cannot swap value as not present in vector");
}

template <typename T>
void VecR<T>::swapValue(T vDel, T vAdd, T vBetween0, T vBetween1) {
    bool swap = false;
    int i;
    int j;
    int k;
    i = 0;
    j = this->n - 1;
    k = 1;
    if (equals(this->v[i], vDel)) {
        if (equals(this->v[j], vBetween0) && equals(this->v[k], vBetween1)) {
            this->v[i] = vAdd;
            swap = true;
        } else if (equals(this->v[j], vBetween1) && equals(this->v[k], vBetween0)) {
            this->v[i] = vAdd;
            swap = true;
        }
    }
    if (swap)
        return;
    for (int i = 1, j = 0, k = 2; i < this->n - 1; ++i, ++j, ++k) {
        if (equals(this->v[i], vDel)) {
            if (equals(this->v[j], vBetween0) && equals(this->v[k], vBetween1)) {
                this->v[i] = vAdd;
                swap = true;
                break;
            } else if (equals(this->v[j], vBetween1) && equals(this->v[k], vBetween0)) {
                this->v[i] = vAdd;
                swap = true;
                break;
            }
        }
    }
    if (swap)
        return;
    i = this->n - 1;
    j = this->n - 2;
    k = 0;
    if (equals(this->v[i], vDel)) {
        if (equals(this->v[j], vBetween0) && equals(this->v[k], vBetween1)) {
            this->v[i] = vAdd;
            swap = true;
        } else if (equals(this->v[j], vBetween1) && equals(this->v[k], vBetween0)) {
            this->v[i] = vAdd;
            swap = true;
        }
    }
    if (swap)
        return;
    std::cout << "Sequence not present in vector" << std::endl;
    throw std::runtime_error("Cannot swap value between values as sequence not present in vector");
}

// Insert value in vector between two others
template <typename T>
void VecR<T>::insertValue(T vInsert, T vBetween0, T vBetween1) {
    if (this->n == this->nMax)
        throw std::runtime_error("Cannot insert to vector to make larger than reserved size");
    bool insert = false;
    int insertPos = -1;

    for (int i = 0; i < this->n; ++i) {
        if (this->v[i] == vBetween0 && this->v[(i + 1) % this->n] == vBetween1) {
            insertPos = (i + 1) % this->n;
            insert = true;
        } else if (this->v[i] == vBetween1 && this->v[(i + 1) % this->n] == vBetween0) {
            insertPos = (i + 1) % this->n;
            insert = true;
        }
        if (insert)
            break;
    }
    if (!insert) {
        for (int i = 0; i < this->n; ++i)
            std::cout << this->v[i] << std::endl;
        throw std::runtime_error("Cannot insert value as surrounding values not present in vector");
    }
    for (int i = this->n; i > insertPos; --i)
        this->v[i] = this->v[i - 1];
    this->v[insertPos] = vInsert;
    ++this->n;
}

// Subscript operator
template <typename T>
T &VecR<T>::operator[](int i) {
    if (this->n <= i)
        throw std::runtime_error("Vector subscript out of bounds");
    return this->v[i];
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
 * @brief Log the vector to the console
 * @param logger The logger to use
 */
template <typename T>
void VecR<T>::toLog(LoggerPtr logger) {
    std::stringstream ss;
    ss << "[";
    for (int i = 0; i < this->n; ++i) {
        ss << this->v[i];
        if (i < this->n - 1)
            ss << ", ";
    }
    ss << "]";
    logger->info(ss.str());
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