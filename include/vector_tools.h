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

template <size_t S>
std::array<double, S> pbcArray(const std::array<double, S> &vector1,
                               const std::array<double, S> &vector2,
                               const std::array<double, S> &dimensions);

template <typename T, size_t S>
  requires Addable<T>
constexpr std::array<T, S> arrayAdd(const std::array<T, S> &array1,
                                    const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Addable<T>
constexpr void arrayAddInPlace(std::array<T, S> &array1,
                               const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> arrayDivide(const std::array<T, S> &array1,
                             const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Divisible<T>
void arrayDivideInPlace(std::array<T, S> &array1,
                        const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> divideArray(const std::array<T, S> &array, const T divisor);

template <typename T, size_t S>
  requires Divisible<T>
void divideArrayInPlace(std::array<T, S> &array, const T divisor);

template <typename T, size_t S>
  requires Divisible<T>
constexpr std::array<T, S> addArray(const std::array<T, S> &array,
                                    const T addition);

template <typename T, size_t S>
  requires Divisible<T>
constexpr void addArray(std::array<T, S> &array, const T addition);

template <typename T>
constexpr std::unordered_set<T>
intersectSets(const std::unordered_set<T> &set1,
              const std::unordered_set<T> &set2);

template <typename Container>
  requires std::ranges::range<Container>
std::string containerToString(const Container &container);

template <typename Container>
  requires Multipliable<typename Container::value_type>
constexpr void
containerMultiply(Container &container,
                  const typename Container::value_type &multiplyBy);

template <typename T>
void setReplace(std::unordered_set<T> &set, const T &oldValue,
                const T &newValue);

template <typename T, size_t S>
  requires Subtractable<T> && Multipliable<T>
T ndDistance(const std::array<T, S> &array1, const std::array<T, S> &array2);

template <typename T, size_t S>
  requires Subtractable<T> && Multipliable<T>
T arrayAbs(const std::array<T, S> &array);

template <typename T, size_t S>
  requires Subtractable<T> && Addable<T>
std::array<T, S> wrapArray(const std::array<T, S> &array,
                           const std::array<T, S> &dimensions);

template <typename T, size_t S>
  requires Subtractable<T> && Addable<T> && Divisible<T>
void wrapArray(std::array<T, S> &array, const std::array<T, S> &dimensions);

template <typename T, size_t S>
  requires Addable<T> && Divisible<T>
std::array<T, S> arrayAverage(const std::vector<std::array<T, S>> &arrays);

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
