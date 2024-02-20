#include "vecf.h"
// ##### VECTOR FIXED SIZE CLASS #####

// Default constructor
template <typename T>
VecF<T>::VecF() {
    this->n = 1;
    this->v = new T[this->n]();
}

// Construct with size
template <typename T>
VecF<T>::VecF(int size) {
    if (size < 0)
        throw std::runtime_error("Vector cannot be instantiated with negative size");
    this->n = size;
    this->v = new T[this->n]();
}

// Construct from VecF
template <typename T>
VecF<T>::VecF(const VecF &source) {
    this->n = source.n;
    this->v = new T[this->n]();
    for (int i = 0; i < this->n; ++i)
        this->v[i] = source.v[i];
}

// Destructor, clear allocated memory
template <typename T>
VecF<T>::~VecF() {
    delete[] this->v;
}

// Check for equality of any value type
template <typename T>
inline bool VecF<T>::equals(const T &a, const T &b) const {
    return a == b || abs(a - b) < abs(std::min(a, b)) * std::numeric_limits<T>::epsilon();
}

// Binary Operators with constant
template <typename T>
void VecF<T>::operator=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] = k;
}

template <typename T>
void VecF<T>::operator+=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] += k;
}

template <typename T>
void VecF<T>::operator-=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] -= k;
}

template <typename T>
void VecF<T>::operator*=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] *= k;
}

template <typename T>
void VecF<T>::operator/=(const T &k) {
    for (int i = 0; i < this->n; ++i)
        this->v[i] /= k;
}

template <typename T>
bool VecF<T>::operator==(const T &k) {
    bool equality = true;
    for (int i = 0; i < this->n; ++i) {
        if (!equals(this->v[i], k))
            equality = false;
    }
    return equality;
}

template <typename T>
bool VecF<T>::operator<(const T &k) {
    bool lt = true;
    for (int i = 0; i < this->n; ++i) {
        if (this->v[i] >= k)
            lt = false;
    }
    return lt;
}

template <typename T>
bool VecF<T>::operator>(const T &k) {
    bool gt = true;
    for (int i = 0; i < this->n; ++i) {
        if (this->v[i] <= k)
            gt = false;
    }
    return gt;
}

template <typename T>
VecF<T> VecF<T>::operator+(const T &k) {
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] + k;
    return vec;
}

template <typename T>
VecF<T> VecF<T>::operator-(const T &k) {
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] - k;
    return vec;
}

template <typename T>
VecF<T> VecF<T>::operator*(const T &k) {
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] * k;
    return vec;
}

template <typename T>
VecF<T> VecF<T>::operator/(const T &k) {
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] / k;
    return vec;
}

// Binary operators with VecF
template <typename T>
void VecF<T>::operator+=(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] += source.v[i];
}

template <typename T>
void VecF<T>::operator-=(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] -= source.v[i];
}

template <typename T>
void VecF<T>::operator*=(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] *= source.v[i];
}

template <typename T>
void VecF<T>::operator/=(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    for (int i = 0; i < this->n; ++i)
        this->v[i] /= source.v[i];
}

template <typename T>
bool VecF<T>::operator==(const VecF &source) {
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
VecF<T> &VecF<T>::operator=(const VecF &source) {
    if (this != &source) {
        delete[] this->v;
        this->n = source.n;
        this->v = new T[this->n]();
        for (int i = 0; i < this->n; ++i)
            this->v[i] = source.v[i];
    }
    return *this;
}

template <typename T>
VecF<T> VecF<T>::operator+(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] + source.v[i];
    return vec;
}

template <typename T>
VecF<T> VecF<T>::operator-(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] - source.v[i];
    return vec;
}

template <typename T>
VecF<T> VecF<T>::operator*(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] * source.v[i];
    return vec;
}

template <typename T>
VecF<T> VecF<T>::operator/(const VecF &source) {
    if (this->n != source.n)
        throw std::runtime_error("Cannot perform binary operation on vectors of different sizes");
    VecF<T> vec(this->n);
    for (int i = 0; i < this->n; ++i)
        vec[i] = this->v[i] / source.v[i];
    return vec;
}

// Unary operators
template <typename T>
VecF<T> VecF<T>::operator-() {
    VecF<T> vec(this->n);
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
bool VecF<T>::contains(T value) {
    for (int i = 0; i < this->n; ++i) {
        if (this->v[i] == value)
            return true;
    }
    return false;
}

template <typename T>
std::string VecF<T>::toString() {
    std::ostringstream os;
    os << *this;
    return os.str();
}

// Subscript operator
template <typename T>
T &VecF<T>::operator[](int i) {
    if (n <= i) {
        std::ostringstream oss;
        oss << "Vector subscript out of bounds: " << i << " >= " << n << " for vector " << *this;
        throw std::runtime_error(oss.str());
    }
    return v[i];
}

// Subscript operator
template <typename T>
T &VecF<T>::operator[](int i) const {
    if (n <= i) {
        std::ostringstream oss;
        oss << "Vector subscript out of bounds: " << i << " >= " << n << " for vector " << *this;
        throw std::runtime_error(oss.str());
    }
    return v[i];
}

/**
 * @brief Stream operator for VecF
 * @param os Output stream
 * @param vec Vector to output
 * @typedef U Type of vector
 * @return Output stream
 */
template <typename U>
std::ostream &operator<<(std::ostream &os, VecF<U> &vec) {
    os << "VecF<" << typeid(U).name() << "> of size " << vec.n << ": [";
    for (const U &value : vec) {
        os << value << ' ';
    }
    os << "]";
    return os;
}

/**
 * @brief Const version of stream operator for VecF
 * @param os Output stream
 * @param vec Vector to output
 * @typedef U Type of vector
 * @return Output stream
 */
template <typename U>
std::ostream &operator<<(std::ostream &os, const VecF<U> &vec) {
    os << "VecF<" << typeid(U).name() << "> of size " << vec.n << ": [";
    for (const U &value : vec) {
        os << value << ' ';
    }
    os << "]";
    return os;
}
