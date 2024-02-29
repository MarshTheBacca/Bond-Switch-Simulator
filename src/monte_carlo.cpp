#include "monte_carlo.h"

// Metropolis-Hastings

// Default constructor
Metropolis::Metropolis() {
    coinFlip = std::uniform_real_distribution<double>(0.0, 1.0);
    randomNumGen.seed(0);
}

// Construct with random seed, temperature and initial energy
Metropolis::Metropolis(int seed, double temperature, double energy) {
    randomNumGen.seed(seed);
    if (temperature <= 0.0)
        throw std::runtime_error("Cannot initialise Metropolis algorithm with zero temperature");
    rTemperature = 1.0 / temperature;
    energyPrev = energy;
}

// Set energy
void Metropolis::setEnergy(double energy) {
    energyPrev = energy;
}

// Get energy
double Metropolis::getEnergy() {
    return energyPrev;
}

// Set temperature
void Metropolis::setTemperature(double temperature) {
    rTemperature = 1.0 / temperature;
}

/**
 * @brief Evaluate Metropolis condition, whether to accept or reject move
 * @param Ef Final energy
 * @param Ei Initial energy
 * @param T_factor Temperature factor
 * @return True if move accepted, false if rejected
 */
bool Metropolis::acceptanceCriterion(double Ef, double Ei, double T_factor) {
    /* Metropolis algorithm efficiently samples Boltzmann distribution
     * 1) if move downhill in energy accept
     * 2) if move uphill accept with probabilitiy min[1,e^-de/t] */
    double deltaE = Ef - Ei;
    if (deltaE < 0.0) {
        return true;
    }
    double probability = exp(-deltaE * rTemperature / T_factor);
    if (coinFlip(randomNumGen) < probability) {
        return true;
    }
    return false;
}
