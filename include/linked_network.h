// Linked reciprocal networks - network and dual pair

#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H
#include "input_data.h"
#include "lammps_object.h"
#include "monte_carlo.h"
#include "network.h"
#include "opt.h"
#include "output_file.h"
#include "pot2d.h"
#include "pot3d.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <random>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string_view>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

enum class SelectionType {
    RANDOM,
    EXPONENTIAL_DECAY
};

class LinkedNetwork {
  public:
    // Data members
    // dual network
    Network networkB;
    int minBCnxs;
    int maxBCnxs;

    // base network
    Network networkA;
    int minACnxs;
    int maxACnxs;

    VecF<double> dimensions;   // Periodic boundary of network, xlo = ylo = 0, so dimensions = [xhi, yhi]
    VecF<double> centreCoords; // Centre of network = [xhi / 2, yhi / 2]

    Network networkT;
    LammpsObject SimpleGraphene;
    LammpsObject TersoffGraphene;
    LammpsObject Triangle_Raft;
    LammpsObject Bilayer;
    LammpsObject BN;

    std::string prefixIn;
    std::string prefixOut;
    bool isOpenMPIEnabled = false;
    bool isSimpleGrapheneEnabled = false;
    bool isTriangleRaftEnabled = false;
    bool isTersoffGrapheneEnabled = false;
    bool isBilayerEnabled = false;
    bool isBNEnabled = false;
    int mcRoutine;
    double CScaling = 1.420;
    double SiScaling = 1.609 * sqrt(32) / 3 / 0.52917721090380;

    VecF<double> crds;       // copy of coordinates in network A (for efficient geometry optimisation)
    std::string mcWeighting; // Either 'weighted' or 'random'
    std::mt19937 mtGen;      // mersenne twister random number generator
    Metropolis mc;           // monte carlo metropolis condition
    VecF<double> potParamsA;
    VecF<double> potParamsB;
    VecF<double> potParamsC;         // potential model parameters (angles, bonds, constraints)
    bool isMaintainConvexityEnabled; // maintain convexity of lattice
    VecF<int> goptParamsA;           // geometry optimisation parameters
    VecF<double> goptParamsB;        // geometry optimisation parameters
    VecF<double> costParams;         // cost function parameters

    VecR<int> fixedRings; // IDs of the fixed rings
    VecR<int> fixedNodes; // IDs of the fixed nodes
    VecR<int> rFixed;

    // Additional data members

    // Constructors
    LinkedNetwork();
    LinkedNetwork(int numRing, LoggerPtr logger);          // Construct hexagonal linked network from scratch
    LinkedNetwork(InputData &inputData, LoggerPtr logger); // Construct from files using an InputData object

    void pushPrefix(std::string_view prefixin, std::string_view prefixout);
    void findFixedRings(bool fixed_rings, std::string filename, LoggerPtr logger);
    void findFixedNodes();

    // Member Functions
    void initialisePotentialModel(double harmonicAngleForceConstant,
                                  double harmonicBondForceConstant,
                                  double harmonicGeometryConstraint,
                                  bool isMaintainConvexityEnabled,
                                  LoggerPtr);                                                  // set up potential model parameters
    void initialiseGeometryOpt(int iterations, double tau, double tolerance, int localExtent); // set up geometry optimsiation parameters
    void initialiseMonteCarlo(const Network &network, double temperature,
                              LoggerPtr logger,
                              int seed = 0);                                         // set up monte carlo
    void initialiseCostFunction(double temperature, int seed, double pk, double rk); // set up cost function
    void rescale(double scaleFactor);                                                // rescale lattice dimensions
    std::tuple<int, int, int, int, int> pickRandomConnection(std::mt19937 &mtGen, SelectionType selectionType);
    int assignValues(int randNodeCoordination, int randNodeConnectionCoordination) const;

    bool generateSwitchIDs(VecF<int> &switchIDsA, VecF<int> &switchIDsB, VecF<int> &switchIDsT,
                           int a, int b, int u, int v, LoggerPtr logger);                                 // get all ids of nodes in lattice A and B needed for switch move
    int findCommonConnection(int idA, int idB, int idDel, LoggerPtr logger);                              //
    int findCommonRing(int idA, int idB, int idDel, LoggerPtr logger);                                    //
    void switchCnx33(VecF<int> switchIDsA, VecF<int> switchIDsB, VecF<int> switchIDsT, LoggerPtr logger); // switch connectivities in lattice
                                                                                                          // between 2x3 coordinate nodes
    bool checkThreeRingEdges(int id);                                                                     // prevent edges being part of three rings
    bool convexRearrangement(VecF<int> switchIDsA, LoggerPtr logger);                                     // rearrange nodes after switch to
                                                                                                          // maintain convexity
    VecF<int> monteCarloSwitchMove(Network network, double &energy, LoggerPtr logger);                    // monte carlo switching move
    VecF<int> monteCarloSwitchMoveLAMMPS(double &SimpleGrapheneEnergy,
                                         double &TersoffGrapheneEnergy,
                                         double &TriangleRaftEnergy,
                                         double &BilayerEnergy, double &BNEnergy,
                                         LoggerPtr logger);
    double globalPotentialEnergy(bool useIntx, bool restrict, Network network, LoggerPtr logger);                             // calculate potential energy of entire system
    VecF<int> globalGeometryOptimisation(bool useIntx, bool restrict, Network network, LoggerPtr logger);                     // geometry optimise entire system
    VecF<int> localGeometryOptimisation(int centreA, int centreB, int extent, bool useIntx, bool restrict, LoggerPtr logger); // geometry optimise subsection of system
    void generateHarmonics(int id, VecR<int> &bonds, VecR<double> &bondParams,
                           VecR<int> &angles, VecR<double> &angleParams,
                           Network network,
                           LoggerPtr logger);                                                            // generate harmonic interactions
    void generateHarmonicsOnly(int id, VecR<int> &bonds, VecR<double> &bondParams, Network network);     // generate harmonic interactions
    void generateRingIntersections(int rId, VecR<int> &intersections);                                   // generate ring intersection interactions
    void generateConvexIntersections(int nId, VecR<int> &intersections);                                 // generate convex ring intersection interactions
    void wrapCoordinates();                                                                              // wrap coordinates if periodic
    void syncCoordinates();                                                                              // update geometry optimised coordinates to networks
    void syncCoordinatesTD();                                                                            // update geometry optimised coordinates to networks
    VecF<double> getNodeDistribution(std::string_view lattice);                                          // get proportion of nodes of each size
    VecF<VecF<int>> getEdgeDistribution(std::string_view lattice);                                       // get unnormalised proportion of node connections
    VecF<double> getAboavWeaire(std::string_view lattice);                                               // get aboav-weaire parameters
    double getAssortativity(std::string_view lattice);                                                   // get network assortativity
    double getAboavWeaireEstimate(std::string_view lattice);                                             // get estimate of aw alpha parameter from assortativity
    VecF<double> getEntropy(std::string_view lattice);                                                   // get node and edge distribution entropy
    VecF<double> getOptimisationGeometry(Network network, VecF<double> &lenHist, VecF<double> &angHist); // get bond/angle mean and standard deviation
    VecF<double> getOptimisationGeometryTD(VecF<double> &lenHist, VecF<double> &angHist);                // get bond/angle mean and standard deviation
    void getRingAreas(VecF<double> &areaSum, VecF<double> &areaSqSum);                                   // get sum of areas and squared
                                                                                                         // areas of each ring size
    double getMaxCluster(std::string_view lattice, int nodeCnd);                                         // get cluster statistics for given node coordination
    VecF<int> getMaxClusters(std::string_view lattice, int minCnd, int maxCnd);                          // get cluster statistics for node coordinations
    bool checkConsistency();                                                                             // check networks are consistent
    bool checkCnxConsistency();                                                                          // check for mutual connections
    bool checkDescriptorConsistency();                                                                   // check descriptors are accurate
    bool checkConvexity();                                                                               // check all angles are convex
    bool checkConvexity(int id);                                                                         // check angles are convex around given node

    // Write Functions
    void writeXYZ(const std::string &prefix);
    void write(const std::string &prefix);
};

#endif // NL_LINKED_NETWORK_H
