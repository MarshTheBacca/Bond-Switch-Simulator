#include "concepts.h"
#include <algorithm>
#include <concepts>
#include <ranges>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

/**
 * @brief Sums all the values in a vector
 * @tparam T The type of the vector
 * @param vector The vector to be summed
 * @return The sum of the vector values
 */
template <typename T> T vectorSum(const std::vector<T> &vector) {
  return std::accumulate(vector.begin(), vector.end(), static_cast<T>(0));
}

/**
 * @brief Converts a vector to a string
 * @tparam T The type of the vector
 * @param vector The vector to be converted
 * @return A string representation of the vector
 */
template <typename T> std::string vectorToString(const std::vector<T> &vector) {
  std::ostringstream oss;
  for (const auto &value : vector) {
    oss << value << " ";
  }
  return oss.str();
}

/**
 * @brief Converts a set to a string
 * @tparam T The type of the set
 * @param set The set to be converted
 * @return A string representation of the set
 */
template <typename T>
std::string setToString(const std::unordered_set<T> &set) {
  std::ostringstream oss;
  for (const auto &value : set) {
    oss << value << " ";
  }
  return oss.str();
}

/**
 * @brief Converts a map to a string
 * @tparam T The type of the keys
 * @tparam U The type of the values
 * @param map The map to be converted
 * @return A string representation of the map
 */
template <typename T, typename U>
std::string mapToString(const std::unordered_map<T, U> &map) {
  std::ostringstream oss;
  for (const auto &value : map) {
    oss << value.first << ": " << value.second << " ";
  }
  return oss.str();
}

/**
 * @brief Calculates the difference between two vectors with periodic boundary
 * conditions
 * @param vector1 First array
 * @param vector2 Second array
 * @param dimensions Dimensions of the system xhi, yhi, (zhi if 3D)
 * @throw std::invalid_argument if the sizes of the arrays do not match the
 * dimensions
 * @return The difference array
 */
template <size_t S>
std::array<double, S> pbcArray(const std::array<double, S> &vector1,
                               const std::array<double, S> &vector2,
                               const std::array<double, S> &dimensions) {
  std::array<double, S> differenceVector;
  double dimensionRange = 0.0;
  double halfDimensionRange = 0.0;
  double difference = 0.0;
  for (size_t i = 0; i < S; ++i) {
    dimensionRange = dimensions[i];
    halfDimensionRange = dimensionRange / 2;
    difference = vector2[i] - vector1[i];
    if (difference > halfDimensionRange) {
      differenceVector[i] = difference - dimensionRange;
    } else if (difference < -halfDimensionRange) {
      differenceVector[i] = difference + dimensionRange;
    }
  }
  return differenceVector;
}

/**
 * @brief Adds two arrays together element-wise
 * @tparam T The type of the arrays
 * @tparam S The size of the arrays
 * @param array1 The first array
 * @param array2 The second array
 * @return An array containing the sum of the two arrays
 */
template <typename T, size_t S>
  requires Addable<T>
std::array<T, S> arrayAdd(const std::array<T, S> &array1,
                          const std::array<T, S> &array2) {
  std::array<T, S> result;
  std::ranges::transform(array1, array2, result.begin(), std::plus<T>());
  return result;
}

/**
 * @brief Subtracts two arrays together element-wise, array2 - array1
 * @tparam T The type of the arrays
 * @tparam S The size of the arrays
 * @param array1 The first array (the subtraction)
 * @param array2 The second array (the base)
 * @return An array containing the difference of the two arrays
 */
template <typename T, size_t S>
  requires Subtractable<T>
std::array<T, S> arraySubtract(const std::array<T, S> &array1,
                               const std::array<T, S> &array2) {
  std::array<T, S> result;
  std::ranges::transform(array2, array1, result.begin(), std::minus<T>());
  return result;
}

/**
 * @brief Adds an array to another array element-wise, IN PLACE
 * @tparam T The type of the arrays
 * @tparam S The size of the arrays
 * @param array1 The first array (WILL BE MODIFIED)
 * @param array2 The second array
 */
template <typename T, size_t S>
  requires Addable<T>
void arrayAdd(std::array<T, S> &array1, const std::array<T, S> &array2) {
  std::ranges::transform(array1, array2, array1.begin(), std::plus<T>());
}

/**
 * @brief Divides an array by another array element-wise
 * @tparam T The type of the arrays
 * @tparam S The size of the arrays
 * @param array1 The first array
 * @param array2 The second array (divisor)
 * @throw std::runtime_error if the divisor contains a zero
 */
template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> arrayDivide(const std::array<T, S> &array1,
                             const std::array<T, S> &array2) {
  std::array<T, S> result;
  std::ranges::transform(array1, array2, result.begin(),
                         [](const T &a, const T &b) {
                           if (b == 0) {
                             throw std::runtime_error("Division by zero");
                           }
                           return a / b;
                         });
  return result;
}

/**
 * @brief Divides an array by another array element-wise, IN PLACE
 * @tparam T The type of the arrays
 * @tparam S The size of the arrays
 * @param array1 The first array (WILL BE MODIFIED)
 * @param array2 The second array (divisor)
 * @throw std::runtime_error if the divisor contains a zero
 */
template <typename T, size_t S>
  requires Divisible<T>
void arrayDivide(std::array<T, S> &array1, const std::array<T, S> &array2) {
  std::ranges::transform(array1, array2, array1.begin(),
                         [](const T &a, const T &b) {
                           if (b == 0) {
                             throw std::runtime_error("Division by zero");
                           }
                           return a / b;
                         });
}

/**
 * @brief Divides an array by a constant
 * @tparam T The type of the array
 * @tparam S The size of the array
 * @param array The array to divide
 * @param divisor The value to divide the array by
 * @throw std::runtime_error if the divisor is zero
 */
template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> divideArray(const std::array<T, S> &array, const T divisor) {
  if (divisor == 0) {
    throw std::runtime_error("Division by zero");
  }
  std::array<T, S> result;
  std::ranges::transform(array, result.begin(), [&divisor](const T &element) {
    return element / divisor;
  });
  return result;
}

/**
 * @brief Divides an array by a constant, IN PLACE
 * @tparam T The type of the array
 * @tparam S The size of the array
 * @param array The array to divide (WILL BE MODIFIED)
 * @param divisor The value to divide the array by
 * @throw std::runtime_error if the divisor is zero
 */
template <typename T, size_t S>
  requires Divisible<T>
void divideArray(std::array<T, S> &array, const T divisor) {
  if (divisor == 0) {
    throw std::runtime_error("Division by zero");
  }
  std::ranges::for_each(array, [&divisor](T &element) { element /= divisor; });
}

/**
 * @brief Adds a constant to every element in an array IN PLACE
 * @tparam T The type of the array
 * @tparam S The size of the array
 * @param array The array to add to
 * @param addition The value to add to the array
 */
template <typename T, size_t S>
  requires Divisible<T>
std::array<T, S> addArray(const std::array<T, S> &array, const T addition) {
  std::array<T, S> result;
  std::ranges::transform(array, result.begin(), [&addition](const T &element) {
    return element + addition;
  });
  return result;
}

/**
 * @brief Adds a constant to every element in an array, IN PLACE
 * @tparam T The type of the array
 * @tparam S The size of the array
 * @param array The array to add to (WILL BE MODIFIED)
 * @param addition The value to add to the array
 */
template <typename T, size_t S>
  requires Divisible<T>
void addArray(std::array<T, S> &array, const T addition) {
  std::ranges::for_each(array,
                        [&addition](T &element) { element += addition; });
}

/**
 * @brief Intersects two sets
 * @tparam T The type of the set elements
 * @param set1 The first set
 * @param set2 The second set
 * @return A set containing the common values
 */
template <typename T>
std::unordered_set<T> intersectSets(const std::unordered_set<T> &set1,
                                    const std::unordered_set<T> &set2) {
  std::unordered_set<T> intersection;
  std::ranges::set_intersection(
      set1, set2, std::inserter(intersection, intersection.begin()));
  return intersection;
}

/**
 * @brief Converts a container to a string
 * @tparam Container The type of the container
 * @param container The container to be converted
 * @return A string representation of the vector
 */
template <typename Container>
  requires std::ranges::range<Container>
std::string containerToString(const Container &container) {
  std::ostringstream oss;
  std::ranges::for_each(container,
                        [&oss](const auto &value) { oss << value << " "; });
  return oss.str();
}

/**
 * @brief Multiplies all the values in a container by a constant
 * @tparam Container The type of the container
 * @param container The container to be multiplied
 * @param multiplyBy The value to multiply the container elements by
 */
template <typename Container>
  requires Multipliable<typename Container::value_type>
void containerMultiply(Container &container,
                       const typename Container::value_type &multiplyBy) {
  std::ranges::for_each(container,
                        [&multiplyBy](auto &value) { value *= multiplyBy; });
}

/**
 * @brief Removes a value and replaces it with another value in a set
 * @tparam T The type of the set elements
 * @param set The set to be modified
 * @param oldValue The value to be replaced
 * @param newValue The value to replace the old value
 */
template <typename T>
void setReplace(std::unordered_set<T> &set, const T &oldValue,
                const T &newValue) {
  set.erase(oldValue);
  set.insert(newValue);
}

/**
 * Calculate the nd distance between two arrays
 * @tparam T The type of the arrays
 * @tparam S The size of the arrays
 * @param array1 The first array
 * @param array2 The second array
 * @return The nd distance between the two arrays
 */
template <typename T, size_t S>
  requires Subtractable<T> && Multipliable<T>
T ndDistance(const std::array<T, S> &array1, const std::array<T, S> &array2) {
  return std::sqrt(std::transform_reduce(
      array1.begin(), array1.end(), array2.begin(), T(0), std::plus<T>(),
      [](T element1, T element2) {
        return (element1 - element2) * (element1 - element2);
      }));
}

template <typename T, size_t S>
  requires Subtractable<T> && Multipliable<T>
T arrayAbs(const std::array<T, S> &array) {
  return std::sqrt(
      std::transform_reduce(array.begin(), array.end(), T(0), std::plus<T>(),
                            [](T element) { return element * element; }));
}
