#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <iostream>
#include <random>

struct Metropolis {
  // Random number generator algorithm
  std::mt19937 randomNumGen = std::mt19937(0);
  // Random number distribution
  std::uniform_real_distribution<double> randNumDist =
      std::uniform_real_distribution<double>(0.0, 1.0);

  /**
   * @brief Default constructor with a random number generator seed of 0
   */
  Metropolis() = default;

  /**
   * @brief Construct with a given seed
   * @param seed seed for random number generator
   */
  explicit Metropolis(const int seed);

  /**
   * @brief Acceptance criterion for Metropolis-Hastings algorithm
   * @param finalEnergy Final energy
   * @param initialEnergy Initial energy
   * @param temperature Temperature factor
   * @return True if final energy less than initial, or with probability
   * e^(-deltaE/T)
   */
  bool acceptanceCriterion(const double finalEnergy, const double initialEnergy,
                           const double temperature);
};

#endif // METROPOLIS_H
