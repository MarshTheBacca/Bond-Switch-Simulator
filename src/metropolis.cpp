#include "metropolis.h"
#include "random_number_generator.h"

bool Metropolis::acceptanceCriterion(const double finalEnergy,
                                     const double initialEnergy,
                                     const double temperature) {
  const double energyChange = finalEnergy - initialEnergy;
  return energyChange < 0 ||
         RandomNumberGenerator::getInstance().getRandomDouble(0.0, 1.0) <
             exp(-energyChange / temperature);
}