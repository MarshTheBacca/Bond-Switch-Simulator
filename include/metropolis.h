#ifndef METROPOLIS_H
#define METROPOLIS_H

#include <iostream>
#include <random>

struct Metropolis {
    std::mt19937 randomNumGen = std::mt19937(0);
    std::uniform_real_distribution<double> randNumDist = std::uniform_real_distribution<double>(0.0, 1.0);

    Metropolis();
    explicit Metropolis(const int &seed);

    bool acceptanceCriterion(const double &finalEnergy, const double &initialEnergy, const double &temperature);
};

#endif // METROPOLIS_H
