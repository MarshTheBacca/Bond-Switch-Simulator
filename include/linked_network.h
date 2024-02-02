// Linked reciprocal networks - network and dual pair

#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H
#include "output_file.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <random>
#include <string_view>
#include "network.h"
#include "pot2d.h"
#include "pot3d.h"
#include "opt.h"
#include "monte_carlo.h"
#include "lammps_object.h"
#include <omp.h>
#include <spdlog/spdlog.h>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

class LinkedNetwork
{
private:
public:
    // Data members
    Network networkA;
    Network networkB; // two reciprocal networks
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
    int MC_Routine;
    double CScaling = 1.420;
    double SiScaling = 1.609 * sqrt((32.0 / 9.0)) / 0.52917721090380;

    VecF<double> crds;       // copy of coordinates in network A (for efficient geometry optimisation)
    std::string MCWeighting; // Either 'weighted' or 'random'
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
    VecR<int> rFixed;
    std::vector<double> weights;

    int spiralRadius; // Radius at which MC moves are considered

    // Additional data members
    int minACnxs;
    int maxACnxs;
    int minBCnxs;
    int maxBCnxs;

    // Constructors
    LinkedNetwork();
    LinkedNetwork(int nodesA, int minCoordination, int maxCoordination, int minRingSize, int maxRingSize, LoggerPtr logger); // construct with starting A lattice
    LinkedNetwork(const std::string &prefixFolderIn, const std::string &prefixFileIn,
                  int minA, int maxA, int minB, int maxB,
                  bool isSimpleGrapheneEnabledArg, bool isTriangleRaftEnabledArg,
                  bool isBilayerEnbledArg, bool isTersoffGrapheneEnabledArg,
                  bool isBNEnabledArg, bool isRestartUsingLAMMPSObjectsEnabledArg,
                  LoggerPtr); // construct by loading from files

    void pushPrefix(std::string_view prefixin, std::string_view prefixout);
    void findFixedRings(bool fixed_rings, std::string filename, LoggerPtr logger);

    void makerFixed();
    // Member Functions
    void initialisePotentialModel(double harmonicAngleForceConstant, double harmonicBondForceConstant,
                                  double harmonicGeometryConstraint, bool isMaintainConvexityEnabled, LoggerPtr); // set up potential model parameters
    void initialiseGeometryOpt(int iterations, double tau, double tolerance, int localExtent);                    // set up geometry optimsiation parameters
    void initialiseMonteCarlo(const Network &network, double temperature, LoggerPtr logger, int seed = 0);        // set up monte carlo
    void initialiseCostFunction(double temperature, int seed, double pk, double rk);                              // set up cost function
    void rescale(double scaleFactor);                                                                             // rescale lattice dimensions
    void project(std::string projType, double param);                                                             // project lattice onto different geometry
    void optimalProjection(std::string projType, LoggerPtr logger);                                               // project lattice onto different geometry with optimal parameters
    int pickSpiralCnx34(int &a, int &b, int &u, int &v, std::mt19937 &gen);
    int pickDiscreteCnx34(int &a, int &b, int &u, int &v, std::mt19937 &gen, LoggerPtr logger);
    int pickRandomCnx34(int &a, int &b, int &u, int &v, std::mt19937 &gen); // choose nodes forming random edge in lattice A, and corresponding nodes in lattice B
    int pickRandomCnx(int &a, int &b, int &u, int &v, std::mt19937 &gen);   // choose nodes forming random edge in lattice A, and corresponding nodes in lattice B
    bool generateSwitchIds34(int cnxType, VecF<int> &switchIdsA, VecF<int> &switchIdsB, VecF<int> &switchIdsT,
                             int a, int b, int u, int v);                                                  // get all ids of nodes in lattice A and B needed for switch move
    int generateMixIds34(int cnxType, VecF<int> &mixIdsA, VecF<int> &mixIdsB, int a, int b, int u, int v); // get all ids of nodes in lattice A and B needed for mix move
    bool generateMixIds(int cnxType, VecF<int> &mixIdsA, VecF<int> &mixIdsB, int a, int b, int u, int v);  // get all ids of nodes in lattice A and B needed for mix move
    int findAssociatedNodeAB(int idA, int idB, int idDel);                                                 //
    int findAssociatedNodeAA(int idA, int idB, int idDel);                                                 //
    void switchCnx33(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);                    // switch connectivities in lattice between 2x3 coordinate nodes
    void switchCnx44(VecF<int> switchIdsA, VecF<int> switchIdsB);                                          // switch connectivities in lattice between 2x4 coordinate nodes
    void switchCnx43(VecF<int> switchIdsA, VecF<int> switchIdsB);                                          // switch connectivities in lattice between 4 and 3 coordinate nodes
    void mixCnx34(VecF<int> mixIdsA, VecF<int> mixIdsB);                                                   // mix connectivities in lattice between 4 and 3 coordinate nodes
    void mixCnx(VecF<int> mixIdsA, VecF<int> mixIdsB);                                                     // mix connectivities in lattice between 4 and 3 coordinate nodes
    bool checkThreeRingEdges(int id);                                                                      // prevent edges being part of three rings
    bool convexRearrangement(int cnxType, VecF<int> switchIdsA, VecF<int> switchIdsB);                     // rearrange nodes after switchi to maintain convexity
    VecF<int> monteCarloSwitchMove(Network network, double &energy, LoggerPtr logger);                     // monte carlo switching move
    VecF<int> SpiralmonteCarloSwitchMoveLAMMPS(int a, double &SimpleGrapheneEnergy,
                                               double &TersoffGrapheneEnergy, double &TriangleRaftEnergy,
                                               double &BilayerEnergy, double &BNEnergy, LoggerPtr logger);
    VecF<int> monteCarloSwitchMoveLAMMPS(double &SimpleGrapheneEnergy, double &TersoffGrapheneEnergy,
                                         double &TriangleRaftEnergy, double &BilayerEnergy,
                                         double &BNEnergy, LoggerPtr logger);
    VecF<int> monteCarloCostSwitchMove(double &cost, double &energy, double pTarget, double rTarget);                         // monte carlo switching move with cost function
    VecF<int> monteCarloMixMove(double &energy, LoggerPtr logger);                                                            // monte carlo mixing move
    double costFunction(double &pTarget, double &rTarget);                                                                    // cost function based on ring statistics and assortative mixing
    double globalPotentialEnergy(bool useIntx, bool restrict, Network network, LoggerPtr logger);                             // calculate potential energy of entire system
    VecF<int> globalGeometryOptimisation(bool useIntx, bool restrict, Network network, LoggerPtr logger);                     // geometry optimise entire system
    VecF<int> localGeometryOptimisation(int centreA, int centreB, int extent, bool useIntx, bool restrict, LoggerPtr logger); // geometry optimise subsection of system
    void generateHarmonics(int id, VecR<int> &bonds, VecR<double> &bondParams,
                           VecR<int> &angles, VecR<double> &angleParams, Network network, LoggerPtr logger); // generate harmonic interactions
    void generateHarmonicsOnly(int id, VecR<int> &bonds, VecR<double> &bondParams, Network network);         // generate harmonic interactions
    void generateRingIntersections(int rId, VecR<int> &intersections);                                       // generate ring intersection interactions
    void generateConvexIntersections(int nId, VecR<int> &intersections);                                     // generate convex ring intersection interactions
    void wrapCoordinates();                                                                                  // wrap coordinates if periodic
    void syncCoordinates();                                                                                  // update geometry optimised coordinates to networks
    void syncCoordinatesTD();                                                                                // update geometry optimised coordinates to networks
    VecF<double> getNodeDistribution(std::string_view lattice);                                              // get proportion of nodes of each size
    VecF<VecF<int>> getEdgeDistribution(std::string_view lattice);                                           // get unnormalised proportion of node connections
    VecF<double> getAboavWeaire(std::string_view lattice);                                                   // get aboav-weaire parameters
    double getAssortativity(std::string_view lattice);                                                       // get network assortativity
    double getAboavWeaireEstimate(std::string_view lattice);                                                 // get estimate of aw alpha parameter from assortativity
    VecF<double> getEntropy(std::string_view lattice);                                                       // get node and edge distribution entropy
    VecF<double> getOptimisationGeometry(Network network, VecF<double> &lenHist, VecF<double> &angHist);     // get bond/angle mean and standard deviation
    VecF<double> getOptimisationGeometryTD(VecF<double> &lenHist, VecF<double> &angHist);                    // get bond/angle mean and standard deviation
    void getRingAreas(VecF<double> &areaSum, VecF<double> &areaSqSum);                                       // get sum of areas and squared areas of each ring size
    double getMaxCluster(std::string_view lattice, int nodeCnd);                                             // get cluster statistics for given node coordination
    VecF<int> getMaxClusters(std::string_view lattice, int minCnd, int maxCnd);                              // get cluster statistics for node coordinations
    bool checkConsistency();                                                                                 // check networks are consistent
    bool checkCnxConsistency();                                                                              // check for mutual connections
    bool checkDescriptorConsistency();                                                                       // check descriptors are accurate
    bool checkConvexity();                                                                                   // check all angles are convex
    bool checkConvexity(int id);                                                                             // check angles are convex around given node

    // Write Functions
    void writeXYZ(const std::string &prefix);
    void write(const std::string &prefix);
};

#endif // NL_LINKED_NETWORK_H
