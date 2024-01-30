// Monte Carlo Methods
#ifndef NL_MONTECARLO_H
#define NL_MONTECARLO_H

#include <iostream>
#include <random>

// Metropolis-Hastings algorithm
class Metropolis
{

private:
    // Data members
    std::mt19937 mtGen;                            // mersenne twister generator
    std::uniform_real_distribution<double> rand01; // uniform distribution
    double rTemperature;                           // reciprocal temperature
    double energyPrev;                             // previous energy

public:
    // Constructors
    Metropolis();
    Metropolis(int seed, double temperature, double energy = 0.0);

    // Member functions
    void setEnergy(double energy);
    void setTemperature(double temperature);
    bool acceptanceCriterion(double Ef, double Ei, double T_factor);
    double getEnergy();
};

#endif // NL_MONTECARLO_H
