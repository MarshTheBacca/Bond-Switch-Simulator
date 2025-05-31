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
 * @brief Gets a vector corresponding to the shortest distance of
 * vector1->vector2 using periodic boundary conditions
 * @param vector1 First array
 * @param vector2 Second array
 * @param dimensions Dimensions of the system xhi, yhi, (zhi if 3D)
 * @note Vectors must be within the period boundary, as no wrapping is applied
 * @return The difference array
 */
template <size_t S>
std::array<double, S> pbcArray(const std::array<double, S> &vector1,
                               const std::array<double, S> &vector2,
                               const std::array<double, S> &dimensions) {
  // Preallocate memory
  std::array<double, S> differenceVector{0.0, 0.0};
  std::array<double, S> halfDimensionRange = divideArray(dimensions, 2.0);
  double difference = 0.0;
  for (size_t i = 0; i < S; ++i) {
    difference = vector2[i] - vector1[i];
    if (difference > halfDimensionRange[i]) {
      differenceVector[i] = difference - dimensions[i];
    } else if (difference < -halfDimensionRange[i]) {
      differenceVector[i] = difference + dimensions[i];
    } else {
      differenceVector[i] = difference;
    }
  }
  return differenceVector;
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
constexpr std::array<T, S> arraySubtract(const std::array<T, S> &array1,
                                         const std::array<T, S> &array2) {
  std::array<T, S> result;
  std::ranges::transform(array2, array1, result.begin(), std::minus<T>());
  return result;
}

/**
 * @brief Adds an array to another array element-wise
 * @tparam T The type of the arrays
 * @tparam S The size of the arrays
 * @param array1 The first array
 * @param array2 The second array
 */
template <typename T, size_t S>
  requires Addable<T>
constexpr std::array<T, S> arrayAdd(const std::array<T, S> &array1,
                                    const std::array<T, S> &array2) {
  std::array<double, 2> result;
  std::ranges::transform(array1, array2, result.begin(), std::plus<T>());
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
constexpr void arrayAddInPlace(std::array<T, S> &array1,
                               const std::array<T, S> &array2) {
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
void arrayDivideInPlace(std::array<T, S> &array1,
                        const std::array<T, S> &array2) {
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
void divideArrayInPlace(std::array<T, S> &array, const T divisor) {
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
constexpr std::array<T, S> addArray(const std::array<T, S> &array,
                                    const T addition) {
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
constexpr void addArray(std::array<T, S> &array, const T addition) {
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
constexpr std::unordered_set<T>
intersectSets(const std::unordered_set<T> &set1,
              const std::unordered_set<T> &set2) {
  std::unordered_set<T> intersection;
  for (const auto &value : set1 | std::views::filter([&set2](const T &v) {
                             return set2.contains(v);
                           })) {
    intersection.insert(value);
  }
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
constexpr void
containerMultiply(Container &container,
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

/**
 * @brief Wraps an array around the periodic boundary conditions
 * @tparam T The type of the array
 * @tparam S The size of the array
 * @param array The array to be wrapped
 * @param dimensions The dimensions of the periodic box
 * @return The wrapped array
 */
template <typename T, size_t S>
  requires Subtractable<T> && Addable<T> && Divisible<T>
std::array<T, S> wrapArray(const std::array<T, S> &array,
                           const std::array<T, S> &dimensions) {
  std::array<T, S> wrappedArray;
  std::transform(array.begin(), array.end(), dimensions.begin(),
                 wrappedArray.begin(), [](const T value, const T dim) {
                   // This is essentially a modulo operation
                   T wrappedValue = value - std::floor(value / dim) * dim;
                   if (wrappedValue < 0) {
                     // If we took off too much, add the dimensions back
                     wrappedValue += dim;
                   }
                   return wrappedValue;
                 });
  return wrappedArray;
}

/**
 * @brief Wraps an array around the periodic boundary conditions IN PLACE
 * @tparam T The type of the array
 * @tparam S The size of the array
 * @param array The array to be wrapped (in place)
 * @param dimensions The dimensions of the periodic box
 * @return The wrapped array
 */
template <typename T, size_t S>
  requires Subtractable<T> && Addable<T> && Divisible<T>
void wrapArray(std::array<T, S> &array, const std::array<T, S> &dimensions) {
  std::transform(array.begin(), array.end(), dimensions.begin(), array.begin(),
                 [](T value, const T dim) {
                   if (dim == 0) {
                     throw std::domain_error("Zero dimension");
                   }
                   // This is essentially a modulo operation
                   T wrappedValue = value - std::floor(value / dim) * dim;
                   if (wrappedValue < 0) {
                     // If we took off too much, add the dimensions back
                     wrappedValue += dim;
                   }
                   return wrappedValue;
                 });
}

template <typename T, size_t S>
  requires Addable<T> && Divisible<T>
std::array<T, S> arrayAverage(const std::vector<std::array<T, S>> &arrays) {
  if (arrays.empty()) {
    throw std::runtime_error("Cannot average an empty vector of arrays");
  }
  std::array<T, S> averageArray = {}; // Initialize to zero
  // Sum all the arrays
  std::ranges::for_each(arrays, [&averageArray](const std::array<T, S> &array) {
    std::ranges::transform(array.begin(), array.end(), averageArray.begin(),
                           averageArray.end(), averageArray.begin(),
                           std::plus<T>());
  });
  // Divide by the number of arrays
  return divideArray(averageArray, static_cast<T>(arrays.size()));
}