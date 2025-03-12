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
 * @brief Checks if a vector contains a value
 * @tparam T The type of the vector
 * @param vector The vector to be checked
 * @param value The value to be checked for
 * @return True if the value is in the vector, false otherwise
 */
template <typename T>
bool vectorContains(const std::vector<T> &vector, const T &value) {
  return std::ranges::find(vector, value) != vector.end();
}
/**
 * @brief Perform linear regression on two vectors
 * @tparam T The type of the vectors
 * @param vector1 The first vector
 * @param vector2 The second vector
 * @return A tuple containing the gradient, intercept and r-squared value
 * @throw std::runtime_error if the vectors are not of the same size
 */
template <typename T>
std::tuple<double, double, double>
vectorLinearRegression(const std::vector<T> &vector1,
                       const std::vector<T> &vector2) {
  if (vector1.size() != vector2.size())
    throw std::runtime_error("Regression error - must be equal number of "
                             "vector1 and vector2 values");

  double sumX = 0.0;
  double sumY = 0.0;
  double sumXY = 0.0;
  double sumXX = 0.0;
  double sumYY = 0.0;

  auto vecSize = static_cast<int>(vector1.size());
  for (int i = 0; i < vecSize; ++i) {
    sumX += vector1[i];
    sumY += vector2[i];
    sumXY += vector1[i] * vector2[i];
    sumXX += vector1[i] * vector1[i];
    sumYY += vector2[i] * vector2[i];
  }

  double sqSumX = sumX * sumX;
  double sqSumY = sumY * sumY;
  double sdX = sqrt(vecSize * sumXX - sqSumX);
  double sdY = sqrt(vecSize * sumYY - sqSumY);
  double r = (vecSize * sumXY - sumX * sumY) / (sdX * sdY);

  // gradient, intercept and r-squared
  double gradient = r * sdY / sdX;
  return std::make_tuple(gradient, (sumY - gradient * sumX) / vecSize, r * r);
}

/**
 * @brief Adds a constant to all the values in a vector
 * @tparam T The type of the vector
 * @param vector The vector to be added to
 * @param addition The value to add to the vector
 */
template <typename T, typename U>
void addToVector(std::vector<T> &vector, const U &addition) {
  for (auto &value : vector) {
    value += addition;
  }
}

/**
 * @brief Subtracts a constant from all the values in a vector
 * @tparam T The type of the vector
 * @tparam U The type of the subtraction value
 * @param vector The vector to be subtracted from
 * @param subtraction The value to subtract from the vector
 */
template <typename T, typename U>
void vectorSubtract(std::vector<T> &vector, const U &subtraction) {
  for (auto &value : vector) {
    value -= static_cast<T>(subtraction);
  }
}

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
 * @brief Multiplies two vectors element-wise
 * @tparam T The type of the vectors
 * @param vector1 The first vector
 * @param vector2 The second vector
 * @return A vector of the multiplied values
 * @throw std::runtime_error if the vectors are not of the same size
 */
template <typename T>
std::vector<T> multiplyVectors(const std::vector<T> &vector1,
                               const std::vector<T> &vector2) {
  auto vecSize = static_cast<int>(vector1.size());
  if (vecSize != vector2.size()) {
    throw std::runtime_error("Vectors must be of the same size");
  }
  std::vector<T> result;
  result.reserve(vecSize);
  for (int i = 0; i < vecSize; ++i) {
    result.push_back(vector1[i] * vector2[i]);
  }
  return result;
}

/**
 * @brief Finds the common values between two vectors
 * @tparam T The type of the vectors
 * @param vector1 The first vector
 * @param vector2 The second vector
 * @return A setunordered_setame T>
 * */
template <typename T>
std::unordered_set<T> intersectVectors(const std::vector<T> &vector1,
                                       const std::vector<T> &vector2) {
  std::unordered_set<T> set2(vector2.begin(), vector2.end());
  std::unordered_set<T> result;
  for (const auto &value : vector1) {
    if (set2.count(value) > 0) {
      result.insert(value);
    }
  }
  return result;
}

/**
 * @brief Deletes all occurrences of a value from a vector
 * @tparam T The type of the vector
 * @param vector The vector from which to delete the value
 * @param args The values to be deleted
 */
template <typename T, typename... Args>
void deleteByValues(std::vector<T> &vector, const Args &...args) {
  std::unordered_set<T> valuesToDelete{args...};
  vector.erase(std::remove_if(vector.begin(), vector.end(),
                              [&valuesToDelete](const T &value) {
                                return valuesToDelete.find(value) !=
                                       valuesToDelete.end();
                              }),
               vector.end());
}

/**
 * @brief Replaces all occurrences of a value in a vector with a new value
 * @tparam T The type of the vector
 * @param vector The vector in which to replace the value
 * @param oldValue The value to be replaced
 * @param newValue The value to replace the old value
 */
template <typename T>
void replaceValue(std::vector<T> &vector, const T &oldValue,
                  const T &newValue) {
  std::ranges::replace(vector, oldValue, newValue);
}

/**
 * @brief Averages the values in a vector
 * @tparam T The type of the vector
 * @param vector The vector to be averaged
 * @return The average of the vector values
 */
template <typename T> double vectorMean(const std::vector<T> &vector) {
  return vectorSum(vector) / vector.size();
}

/**
 * @brief Gets the unique values in a vector
 * @tparam T The type of the vector
 * @param vector The vector to be checked
 * @return A vector of unique values
 */
template <typename T>
std::vector<T> getUniqueValues(const std::vector<T> &vector) {
  std::unordered_set<T> uniqueSet;
  std::vector<T> uniqueVector;
  for (const auto &value : vector) {
    if (uniqueSet.insert(value).second) {
      uniqueVector.push_back(value);
    }
  }
  return uniqueVector;
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
 * @brief Adds values to a vector
 * @tparam T The type of the vector
 * @tparam Args The type of the values to be added
 * @param vector The vector to which the values will be added
 */
template <typename T, typename... Args>
void vectorAddValues(std::vector<T> &vector, const Args &...args) {
  (vector.push_back(args), ...);
}

/**
 * @brief Combines two vectors into one
 * @tparam T The type of the vectors
 * @param vec1 The first vector
 * @param vec2 The second vector
 * @return A vector containing all the values from both vectors
 */
template <typename T>
std::vector<T> combineVectors(const std::vector<T> &vec1,
                              const std::vector<T> &vec2) {
  std::vector<T> result;
  result.reserve(vec1.size() +
                 vec2.size()); // Reserve enough space for all elements
  result.insert(result.end(), vec1.begin(), vec1.end());
  result.insert(result.end(), vec2.begin(), vec2.end());
  return result;
}

/**
 * @brief Normalises a vector in place
 * @tparam T The type of the vector
 * @param vec The vector to be normalised
 */
template <typename T> void vectorNormalise(std::vector<T> &vec) {
  vectorDivide(vec, vectorSum(vec));
}

/**
 * @brief Finds the difference between two sets
 * @tparam T The type of the sets
 * @param set1 The first set
 * @param set2 The second set
 * @return A set containing the values in set1 that are not in set2
 */
template <typename T>
std::unordered_set<T> setDifference(const std::unordered_set<T> &set1,
                                    const std::unordered_set<T> &set2) {
  std::unordered_set<T> result;
  for (const auto &value : set1) {
    if (set2.count(value) == 0) {
      result.insert(value);
    }
  }
  return result;
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
 * @brief Gets the value at a specific index in a set
 * @tparam T The type of the set elements
 * @param set The set to be accessed
 * @param index The index of the value to be accessed
 */
template <typename T> T &getSetAt(std::set<T> &set, const size_t index) {
  auto it = set.begin();
  std::advance(it, index);
  return *it;
}
