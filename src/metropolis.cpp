#include "metropolis.h"

/**
 * @brief Default constructor with a random number generator seed of 0
 */
Metropolis::Metropolis() = default;

/**
 * @brief Construct with a given seed
 * @param seed seed for random number generator
 */
Metropolis::Metropolis(const int &seed) {
    randomNumGen.seed(seed);
}

/**
 * @brief Acceptance criterion for Metropolis-Hastings algorithm
 * @param finalEnergy Final energy
 * @param initialEnergy Initial energy
 * @param temperature Temperature factor
 * @return True if final energy less than initial, or with probability e^(-deltaE/T)
 */
bool Metropolis::acceptanceCriterion(const double &finalEnergy, const double &initialEnergy, const double &temperature) {
    const double energyChange = finalEnergy - initialEnergy;
    return energyChange < 0 || randNumDist(randomNumGen) < exp(-energyChange / temperature);
}