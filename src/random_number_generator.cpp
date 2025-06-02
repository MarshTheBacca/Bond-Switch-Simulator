#include "random_number_generator.h"
#include <random>

int RandomNumberGenerator::getRandomInt(const int min, const int max) {
  return this->getInstance()._getRandomInt(min, max);
}

int RandomNumberGenerator::_getRandomInt(const int min, const int max) {
  this->intDist.param(std::uniform_int_distribution<int>::param_type(min, max));
  return this->intDist(rng);
}

double RandomNumberGenerator::getRandomDouble(const double min,
                                              const double max) {
  return this->getInstance()._getRandomDouble(min, max);
}

double RandomNumberGenerator::_getRandomDouble(const double min,
                                               const double max) {
  this->doubleDist.param(
      std::uniform_real_distribution<double>::param_type(min, max));
  return this->doubleDist(rng);
}
