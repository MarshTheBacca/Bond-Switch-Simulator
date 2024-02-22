// Linked reciprocal networks - network and dual pair

#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H
#include "input_data.h"
#include "lammps_object.h"
#include "monte_carlo.h"
#include "network.h"
#include "opt.h"
#include "output_file.h"
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

class LinkedNetwork {
  public:
    // Data members

    Network networkB; // Ring network
    int minBCnxs;     // Minimum coordination number of ring network
    int maxBCnxs;     // Maximum coordination number of ring network

    Network networkA; // Base network
    int minACnxs;     // Minimum coordination number of base network
    int maxACnxs;     // Maximum coordination number of base network

    VecF<double> dimensions;   // Periodic boundary of network, xlo = ylo = 0, so dimensions = [xhi, yhi]
    VecF<double> centreCoords; // Centre of network = [xhi / 2, yhi / 2]

    LammpsObject lammpsNetwork;  // LAMMPS object for network
    double energy;               // The current energy of the system
    int numSwitches = 0;         // Number of switches performed
    int numAcceptedSwitches = 0; // Number of switches accepted
    std::vector<double> currentCoords;

    std::string prefixIn;
    std::string prefixOut;

    bool isOpenMPIEnabled;       // Whether to use MPI
    SelectionType selectionType; // Either 'weighted' or 'random'
    std::mt19937 mtGen;          // mersenne twister random number generator
    Metropolis mc;               // monte carlo metropolis condition
    double weightedDecay = 1;    // decay factor for weighted monte carlo
    double maximumBondLength;    // Maximum bond length
    double maximumAngle;         // Maximum angle between atoms
    bool writeMovie;             // Write movie file

    std::vector<int> fixedRings; // IDs of the fixed rings
    std::vector<int> fixedNodes; // IDs of the fixed nodes

    int failedBondLengthChecks = 0; // Number of failed bond length checks
    int failedAngleChecks = 0;      // Number of failed angle checks
    int failedEnergyChecks = 0;     // Number of failed energy checks

    // Additional data members

    // Constructors
    LinkedNetwork();
    LinkedNetwork(int numRing, LoggerPtr logger);                // Construct hexagonal linked network from scratch
    LinkedNetwork(const InputData &inputData, LoggerPtr logger); // Construct from files using an InputData object

    void findFixedRings(const std::string &filename, LoggerPtr logger);
    void findFixedNodes();

    // Member Functions
    void rescale(double scaleFactor); // rescale lattice dimensions
    std::tuple<int, int, int, int> pickRandomConnection(std::mt19937 &mtGen, const SelectionType &selectionType);
    int assignValues(int randNodeCoordination, int randNodeConnectionCoordination) const;

    int findCommonConnection(int idA, int idB, int idDel, LoggerPtr logger); //
    int findCommonRing(int idA, int idB, int idDel, LoggerPtr logger);       //

    // between 2x3 coordinate nodes
    bool checkThreeRingEdges(const int &rinNode) const;               // prevent edges being part of three rings
    bool convexRearrangement(VecF<int> switchIDsA, LoggerPtr logger); // rearrange nodes after switch to
                                                                      // maintain convexity
    void monteCarloSwitchMoveLAMMPS(LoggerPtr logger);
    void generateHarmonics(int id, VecR<int> &bonds, VecR<double> &bondParams,
                           VecR<int> &angles, VecR<double> &angleParams,
                           Network network,
                           LoggerPtr logger);                                                        // generate harmonic interactions
    void generateHarmonicsOnly(int id, VecR<int> &bonds, VecR<double> &bondParams, Network network); // generate harmonic interactions
    void generateRingIntersections(int rId, VecR<int> &intersections);                               // generate ring intersection interactions
    void generateConvexIntersections(int nId, VecR<int> &intersections);
    void wrapCoords(std::vector<double> &coords);                               // generate convex ring intersection interactions
    VecF<double> getNodeDistribution(std::string_view lattice);                 // get proportion of nodes of each size
    VecF<VecF<int>> getEdgeDistribution(std::string_view lattice);              // get unnormalised proportion of node connections
    VecF<double> getAboavWeaire(std::string_view lattice);                      // get aboav-weaire parameters
    double getAssortativity(std::string_view lattice);                          // get network assortativity
    double getAboavWeaireEstimate(std::string_view lattice);                    // get estimate of aw alpha parameter from assortativity
    VecF<double> getEntropy(std::string_view lattice);                          // get node and edge distribution entropy
    void getRingAreas(VecF<double> &areaSum, VecF<double> &areaSqSum);          // get sum of areas and squared
                                                                                // areas of each ring size
    double getMaxCluster(std::string_view lattice, int nodeCnd);                // get cluster statistics for given node coordination
    VecF<int> getMaxClusters(std::string_view lattice, int minCnd, int maxCnd); // get cluster statistics for node coordinations
    bool checkConsistency();                                                    // check networks are consistent
    bool checkCnxConsistency();                                                 // check for mutual connections
    bool checkDescriptorConsistency();                                          // check descriptors are accurate

    // Write Functions
    void writeXYZ(const std::string &prefix);
    void write(const std::string &prefix);

    void pushCoords(std::vector<double> &coords);
    void showCoords(std::vector<double> &coords, LoggerPtr logger) const;

    bool genSwitchOperations(int baseNode1, int baseNode2, int ringNode1, int ringNode2,
                             std::vector<int> &bondBreaks, std::vector<int> &bondMakes,
                             std::vector<int> &angleBreaks, std::vector<int> &angleMakes,
                             std::vector<int> &ringBondBreakMake, std::vector<int> &convexCheckIDs, LoggerPtr logger);

    void switchNetMCGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &ringBondBreakMake);
    void revertNetMCGraphene(const std::vector<Node> &initialInvolvedNodesA, const std::vector<Node> &initialInvolvedNodesB);

    std::tuple<std::vector<double>, std::vector<double>> rotateBond(const int &atomID1, const int &atomID2,
                                                                    const Direction &direct, LoggerPtr logger) const;
    Direction getRingsDirection(const std::vector<int> &ringNodeIDs, LoggerPtr logger) const;

    bool checkClockwiseNeighbours(const int &nodeID) const;
    bool checkClockwiseNeighbours(const int &nodeID, const std::vector<double> &coords) const;

    bool checkAllClockwiseNeighbours(LoggerPtr logger) const;

    bool checkAnglesWithinRange(const std::vector<double> &coords, LoggerPtr logger);
    bool checkAnglesWithinRange(const std::vector<int> &nodeIDs, const std::vector<double> &coords, LoggerPtr logger);
    void arrangeNeighboursClockwise(const int &nodeID, const std::vector<double> &coords);
    bool checkBondLengths(const int &nodeID, const std::vector<double> &coords, LoggerPtr logger) const;
    bool checkBondLengths(const std::vector<int> &nodeIDs, const std::vector<double> &coords, LoggerPtr logger) const;
};

#endif // NL_LINKED_NETWORK_H
