#include "metropolis.h"

Metropolis::Metropolis(const int seed) : randomNumGen(seed) {}

bool Metropolis::acceptanceCriterion(const double finalEnergy,
                                     const double initialEnergy,
                                     const double temperature) {
  const double energyChange = finalEnergy - initialEnergy;
  return energyChange < 0 ||
         randNumDist(randomNumGen) < exp(-energyChange / temperature);
}