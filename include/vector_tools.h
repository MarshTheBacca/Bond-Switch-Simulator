// Functions to Operate on Vector Classes

#ifndef VEC_TOOLS_H
#define VEC_TOOLS_H

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>
#include <unordered_set>
#include <vector>

// template functions

template <typename T>
bool vectorContains(const std::vector<T> &vector, const T &value);

template <typename T>
std::tuple<double, double, double> vectorLinearRegression(const std::vector<T> &vector1, const std::vector<T> &vector2);

template <typename T>
void vectorDivide(std::vector<T> &vector, const double &divideBy);

template <typename T>
void vectorMultiply(std::vector<T> &vector, const double &multiplyBy);

template <typename T, typename U>
void addToVector(std::vector<double> &vector, const U &addition);

template <typename T, typename U>
void vectorSubtract(std::vector<T> &vector, const U &subtraction);

template <typename T>
T vectorSum(const std::vector<T> &vector);

template <typename T>
std::vector<T> multiplyVectors(const std::vector<T> &vector1, const std::vector<T> &vector2);

template <typename T>
std::unordered_set<T> intersectVectors(const std::vector<T> &vector1, const std::vector<T> &vector2);

template <typename T, typename... Args>
void deleteByValues(std::vector<T> &vector, const Args &...args);

template <typename T>
void replaceValue(std::vector<T> &vector, const T &oldValue, const T &newValue);

template <typename T>
double vectorMean(const std::vector<T> &vector);

template <typename T>
std::vector<T> getUniqueValues(const std::vector<T> &vector);

template <typename T>
std::string vectorToString(const std::vector<T> &vector);

template <typename T>
std::string setToString(const std::unordered_set<T> &set);

template <typename T, typename... Args>
void vectorAddValues(std::vector<T> &vector, const Args &...args);

template <typename T>
std::vector<T> combineVectors(const std::vector<T> &vec1, const std::vector<T> &vec2);

template <typename T>
void vectorNormalise(const std::vector<T> &vec);

template <typename T>
std::unordered_set<T> setDifference(const std::unordered_set<T> &set1, const std::unordered_set<T> &set2);

// Regular functions

std::vector<double> pbcVector(const std::vector<double> &vector1, const std::vector<double> &vector2, const std::vector<double> &dimensions);
void normaliseMap(std::map<int, double> &map);
void showNestedMap(const std::map<int, std::map<int, double>> &map);
double getClockwiseAngleBetweenVectors(const std::vector<double> &vector1, const std::vector<double> &vector2);
double getClockwiseAngle(const std::vector<double> &coord1, const std::vector<double> &coord2, const std::vector<double> &dimensions);
double getClockwiseAngle(const std::vector<double> &point);

void sortCoordinatesClockwise(std::vector<std::vector<double>> &coords);
double calculatePolygonArea(std::vector<std::vector<double>> &vertices);

#include "vector_tools.tpp"
#endif // VEC_TOOLS_H
