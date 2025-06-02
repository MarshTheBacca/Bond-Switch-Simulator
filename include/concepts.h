
#ifndef CONCEPTS_H
#define CONCEPTS_H
#include <concepts>

// Concept to check if a type supports the multiply operator
template <typename T>
concept Multipliable = requires(T a, const T b) {
  { a *= b } -> std::same_as<T &>;
};

// Concept to check if a type supports the divide operator
template <typename T>
concept Divisible = requires(T a, const T b) {
  { a /= b } -> std::same_as<T &>;
};

// Concept to check if a type supports the addition operator
template <typename T>
concept Addable = requires(T a, const T b) {
  { a += b } -> std::same_as<T &>;
};

// Concept to check if a type supports the subtraction operator
template <typename T>
concept Subtractable = requires(T a, const T b) {
  { a -= b } -> std::same_as<T &>;
};

#endif // CONCEPTS_H