#ifndef VEC_TOOLS_H
#define VEC_TOOLS_H

#include "concepts.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <unordered_set>
#include <vector>

// template functions

template <typename T>
bool vectorContains(const std::vector<T> &vector, const T &value);

template <typename T>
std::tuple<double, double, double>
vectorLinearRegression(const std::vector<T> &vector1,
                       const std::vector<T> &vector2);

template <typename T>
void vectorDivide(std::vector<T> &vector, const double &divideBy);

template <typename T, typename U>
void addToVector(std::vector<double> &vector, const U &addition);

template <typename T, typename U>
void vectorSubtract(std::vector<T> &vector, const U &subtraction);

template <typename T> T vectorSum(const std::vector<T> &vector);

template <typename T>
std::vector<T> multiplyVectors(const std::vector<T> &vector1,
                               const std::vector<T> &vector2);

template <typename T>
std::unordered_set<T> intersectVectors(const std::vector<T> &vector1,
                                       const std::vector<T> &vector2);

template <typename T, typename... Args>
void deleteByValues(std::vector<T> &vector, const Args &...args);

template <typename T>
void replaceValue(std::vector<T> &vector, const T &oldValue, const T &newValue);

template <typename T> double vectorMean(const std::vector<T> &vector);

template <typename T>
std::vector<T> getUniqueValues(const std::vector<T> &vector);

template <typename T> std::string vectorToString(const std::vector<T> &vector);

template <typename T> std::string setToString(const std::unordered_set<T> &set);

template <typename T, typename... Args>
void vectorAddValues(std::vector<T> &vector, const Args &...args);

template <typename T>
std::vector<T> combineVectors(const std::vector<T> &vec1,
                              const std::vector<T> &vec2);

template <typename T> void vectorNormalise(const std::vector<T> &vec);

template <typename T>
std::unordered_set<T> setDifference(const std::unordered_set<T> &set1,
                                    const std::unordered_set<T> &set2);

template <size_t S>
std::array<double, S> pbcArray(const std::array<double, S> &vector1,
                               const std::array<double, S> &vector2,
                               const std::array<double, S> &dimensions);

template <typename T, size_t S>
  requires Addable<T>
std::array<T, S> arrayAdd(const std::array<T, S> &array1,
                          const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Addable<T>
void arrayAdd(std::array<T, S> &array1, const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> arrayDivide(const std::array<T, S> &array1,
                             const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Divisible<T>
void arrayDivide(std::array<T, S> &array1, const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> divideArray(const std::array<T, S> &array, const T divisor);

template <typename T, size_t S>
  requires Divisible<T>
void divideArray(std::array<T, S> &array, const T divisor);

template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> addArray(const std::array<T, S> &array, const T addition);

template <typename T, size_t S>
  requires Divisible<T>
void addArray(std::array<T, S> &array, const T addition);

template <typename T>
std::set<T> intersectSets(const std::set<T> &set1, const std::set<T> &set2);

template <typename Container>
  requires std::ranges::range<Container>
std::string containerToString(const Container &container);

template <typename Container>
  requires Multipliable<typename Container::value_type>
void containerMultiply(Container &container,
                       const typename Container::value_type &multiplyBy);

template <typename T>
void setReplace(std::set<T> &set, const T &oldValue, const T &newValue);

template <typename T> T &getSetAt(std::set<T> &set, size_t index);

// Regular functions
void normaliseMap(std::map<int, double> &map);
void showNestedMap(const std::map<int, std::map<int, double>> &map);
double getClockwiseAngleBetweenVectors(const std::array<double, 2> &vector1,
                                       const std::array<double, 2> &vector2);
double getClockwiseAngle(const std::array<double, 2> &coord1,
                         const std::array<double, 2> &coord2,
                         const std::array<double, 2> &dimensions);
double getClockwiseAngle(const std::array<double, 2> &point);

void sortCoordinatesClockwise(std::vector<std::array<double, 2>> &coords);
double calculatePolygonArea(std::vector<std::array<double, 2>> &vertices);

#include "vector_tools.tpp"
#endif // VEC_TOOLS_H
