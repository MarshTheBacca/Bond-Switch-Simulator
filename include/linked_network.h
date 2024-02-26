// Linked reciprocal networks - network and dual pair

#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H
#include "input_data.h"
#include "lammps_object.h"
#include "monte_carlo.h"
#include "network.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string_view>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

enum class Direction {
    CLOCKWISE,
    ANTICLOCKWISE
};

struct LinkedNetwork {
    // Data members

    Network networkB; // Ring network
    int minBCnxs;     // Minimum coordination number of ring network
    int maxBCnxs;     // Maximum coordination number of ring network

    Network networkA; // Base network
    int minACnxs;     // Minimum coordination number of base network
    int maxACnxs;     // Maximum coordination number of base network
    int analysisMinBCnxs;
    int analysisMaxBCnxs;

    std::vector<double> dimensions;   // Periodic boundary of network, xlo = ylo = 0, so dimensions = [xhi, yhi]
    std::vector<double> centreCoords; // Centre of network = [xhi / 2, yhi / 2]

    LammpsObject lammpsNetwork; // LAMMPS object for network
    double energy;              // The current energy of the system

    std::vector<double> currentCoords;

    bool isOpenMPIEnabled;          // Whether to use MPI
    SelectionType selectionType;    // Either 'weighted' or 'random'
    std::mt19937 mtGen;             // mersenne twister random number generator
    Metropolis metropolisCondition; // monte carlo metropolis condition
    double weightedDecay;           // decay factor for weighted monte carlo
    double maximumBondLength;       // Maximum bond length
    double maximumAngle;            // Maximum angle between atoms
    bool writeMovie;                // Write movie file or not

    std::unordered_set<int> fixedRings; // IDs of the fixed rings
    std::unordered_set<int> fixedNodes; // IDs of the fixed nodes

    int numSwitches = 0;            // Number of switches performed
    int numAcceptedSwitches = 0;    // Number of switches accepted
    int failedBondLengthChecks = 0; // Number of failed bond length checks
    int failedAngleChecks = 0;      // Number of failed angle checks
    int failedEnergyChecks = 0;     // Number of failed energy checks

    LoggerPtr logger; // Logger
    std::vector<double> weights;

    // Constructors
    LinkedNetwork();
    LinkedNetwork(const int &numRing, LoggerPtr logger);         // Construct hexagonal linked network from scratch
    LinkedNetwork(const InputData &inputData, LoggerPtr logger); // Construct from files using an InputData object

    void findFixedRings(const std::string &filename);
    void findFixedNodes();

    // Member Functions
    void rescale(double scaleFactor); // rescale lattice dimensions
    void updateWeights();
    std::tuple<int, int, int, int> pickRandomConnection();
    int assignValues(int randNodeCoordination, int randNodeConnectionCoordination) const;

    int findCommonConnection(const int &baseNode, const int &ringNode, const int &excludeNode) const;
    int findCommonRing(const int &baseNode1, const int &baseNode2, const int &excludeNode) const;

    void monteCarloSwitchMoveLAMMPS();
    bool checkConsistency();           // check networks are consistent
    bool checkCnxConsistency();        // check for mutual connections
    bool checkDescriptorConsistency(); // check descriptors are accurate

    void write(const std::string &prefix);

    void pushCoords(const std::vector<double> &coords);
    void showCoords(const std::vector<double> &coords) const;
    void wrapCoords(std::vector<double> &coords) const; // generate convex ring intersection interactions

    bool genSwitchOperations(int baseNode1, int baseNode2, int ringNode1, int ringNode2,
                             std::vector<int> &bondBreaks, std::vector<int> &bondMakes,
                             std::vector<int> &angleBreaks, std::vector<int> &angleMakes,
                             std::vector<int> &ringBondBreakMake, std::vector<int> &convexCheckIDs);

    void switchNetMCGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &ringBondBreakMake);
    void revertNetMCGraphene(const std::vector<Node> &initialInvolvedNodesA, const std::vector<Node> &initialInvolvedNodesB);

    std::tuple<std::vector<double>, std::vector<double>> rotateBond(const int &atomID1, const int &atomID2,
                                                                    const Direction &direct) const;
    Direction getRingsDirection(const std::vector<int> &ringNodeIDs) const;

    bool checkClockwiseNeighbours(const int &nodeID) const;
    bool checkClockwiseNeighbours(const int &nodeID, const std::vector<double> &coords) const;
    bool checkAllClockwiseNeighbours() const;
    void arrangeNeighboursClockwise(const int &nodeID, const std::vector<double> &coords);

    bool checkAnglesWithinRange(const std::vector<double> &coords);
    bool checkAnglesWithinRange(const std::vector<int> &nodeIDs, const std::vector<double> &coords);
    bool checkBondLengths(const int &nodeID, const std::vector<double> &coords) const;
    bool checkBondLengths(const std::vector<int> &nodeIDs, const std::vector<double> &coords) const;

    std::vector<double> getRingSizes() const;
    double calculatePolygonArea(const std::vector<std::vector<double>> &vertices) const;
    std::vector<double> getRingAreas() const;
};

#endif // NL_LINKED_NETWORK_H
