#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <iostream>
#include <random>

struct Metropolis {

  /**
   * @brief Default constructor with a random number generator seed of 0
   */
  Metropolis() = default;

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
