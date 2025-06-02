#ifndef RANDOM_NUMBER_GENERATOR_H
#define RANDOM_NUMBER_GENERATOR_H

#include <concepts>
#include <random>

class RandomNumberGenerator {
public:
  // Delete copy constructor and assignment operator
  RandomNumberGenerator(const RandomNumberGenerator &) = delete;
  RandomNumberGenerator &operator=(const RandomNumberGenerator &) = delete;

  // Static method to initialize the singleton instance with a seed
  static void initialize(const unsigned int seed) {
    getInstance().rng.seed(seed);
  }

  // Static method to get the singleton instance
  static RandomNumberGenerator &getInstance() {
    static RandomNumberGenerator instance;
    return instance;
  }

  int getRandomInt(const int min, const int max);

  double getRandomDouble(const double min, const double max);

  template <std::input_iterator T>
  T getRandomElement(const T begin, const T end) {
    return getInstance()._getRandomElement(begin, end);
  }

private:
  std::mt19937 rng;
  std::uniform_int_distribution<int> intDist;
  std::uniform_real_distribution<double> doubleDist;
  std::uniform_int_distribution<size_t> sizeDist;

  // Private constructor
  RandomNumberGenerator() : rng(std::random_device{}()) {}

  int _getRandomInt(const int min, const int max);
  double _getRandomDouble(const double min, const double max);

  template <std::input_iterator T>
  T _getRandomElement(const T begin, const T end) {
    size_t range = std::distance(begin, end) - 1;
    sizeDist.param(std::uniform_int_distribution<size_t>::param_type(0, range));
    return std::next(begin, sizeDist(rng));
  }
};

#endif // RANDOM_NUMBER_GENERATOR_H