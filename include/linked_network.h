// Linked reciprocal networks - network and dual pair

#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H
#include "output_file.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <random>
#include "network.h"
#include "pot2d.h"
#include "pot3d.h"
#include "opt.h"
#include "monte_carlo.h"
#include "lammps_c_interface.h"
#include <omp.h>
#include <spdlog/spdlog.h>

using namespace std;
using LoggerPtr = std::shared_ptr<spdlog::logger>;

class LinkedNetwork
{
private:
public:
    // Data members
    Network networkA, networkB; // two reciprocal networks
    Network networkT;
    LammpsObject SimpleGraphene;
    LammpsObject TersoffGraphene;
    LammpsObject Triangle_Raft;
    LammpsObject Bilayer;
    LammpsObject BN;

    string prefixIn, prefixOut;
    bool isOpenMP = false;
    bool isSimpleGraphene = true;
    bool isTriangleRaft = false;
    bool isTersoffGraphene = false;
    bool isBilayer = false;
    bool isBN = false;
    int MC_Routine;
    double CScaling = 1.420;
    double SiScaling = 1.609 * sqrt((32.0 / 9.0)) / 0.52917721090380;
    //    LammpsObject Bilayer;

    VecF<double> crds; // copy of coordinates in network A (for efficient geometry optimisation)
    string MCWeighting;
    mt19937 mtGen;                                   // mersenne twister random number generator
    Metropolis mc, mcCost;                           // monte carlo metropolis condition
    VecF<double> potParamsA, potParamsB, potParamsC; // potential model parameters (angles, bonds, constraints)
    VecF<int> potParamsD;                            // potential model intersection/convex parameters
    VecF<int> goptParamsA;                           // geometry optimisation parameters
    VecF<double> goptParamsB;                        // geometry optimisation parameters
    VecF<double> costParams;                         // cost function parameters

    VecR<int> fixedRings;
    VecR<int> rFixed;
    vector<double> weights;

    int spiralRadius;

    // Additional data members
    int minACnxs, maxACnxs, minBCnxs, maxBCnxs;

    // Constructors
    LinkedNetwork();
    LinkedNetwork(int nodesA, int minA, int maxA, int minB, int maxB); // construct with starting A lattice
    LinkedNetwork(string prefixFolderIn, string prefixFileIn, string prefixFolderOut,
                  int minA, int maxA, int minB, int maxB,
                  bool isSimpleGraphene, bool isTriangleRaft, bool isBilayer, bool isTersoffGraphene, bool isBN, bool restartLammps, LoggerPtr); // construct by loading from files

    void pushPrefix(string prefixin, string prefixout);
    void findFixedRings(bool fixed_rings, string filename, LoggerPtr logger);

    void makerFixed();
    // Member Functions
    void initialisePotentialModel(Network network,
                                  double harmonicAngleForceConstant, double harmonicBondForceConstant,
                                  double harmonicGeometryConstraint, bool isMaintainConvexityEnabled, LoggerPtr); // set up potential model parameters
    void initialiseGeometryOpt(int iterations, double tau, double tolerance, int localExtent);                    // set up geometry optimsiation parameters
    void initialiseMonteCarlo(Network network, double temperature, int seed = 0);                                 // set up monte carlo
    void initialiseCostFunction(double temperature, int seed, double pk, double rk);                              // set up cost function
    void rescale(double scaleFactor);                                                                             // rescale lattice dimensions
    void project(string projType, double param);                                                                  // project lattice onto different geometry
    void optimalProjection(string projType);                                                                      // project lattice onto different geometry with optimal parameters
    int pickSpiralCnx34(int &a, int &b, int &u, int &v, mt19937 &gen);
    int pickDiscreteCnx34(int &a, int &b, int &u, int &v, mt19937 &gen);
    int pickRandomCnx34(int &a, int &b, int &u, int &v, mt19937 &gen);                                                                     // choose nodes forming random edge in lattice A, and corresponding nodes in lattice B
    int pickRandomCnx(int &a, int &b, int &u, int &v, mt19937 &gen);                                                                       // choose nodes forming random edge in lattice A, and corresponding nodes in lattice B
    int generateSwitchIds34(int cnxType, VecF<int> &switchIdsA, VecF<int> &switchIdsB, VecF<int> &switchIdsT, int a, int b, int u, int v); // get all ids of nodes in lattice A and B needed for switch move
    int generateMixIds34(int cnxType, VecF<int> &mixIdsA, VecF<int> &mixIdsB, int a, int b, int u, int v);                                 // get all ids of nodes in lattice A and B needed for mix move
    int generateMixIds(int cnxType, VecF<int> &mixIdsA, VecF<int> &mixIdsB, int a, int b, int u, int v);                                   // get all ids of nodes in lattice A and B needed for mix move
    int findAssociatedNodeAB(int idA, int idB, int idDel);                                                                                 //
    int findAssociatedNodeAA(int idA, int idB, int idDel);                                                                                 //
    void switchCnx33(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);                                                    // switch connectivities in lattice between 2x3 coordinate nodes
    void switchCnx44(VecF<int> switchIdsA, VecF<int> switchIdsB);                                                                          // switch connectivities in lattice between 2x4 coordinate nodes
    void switchCnx43(VecF<int> switchIdsA, VecF<int> switchIdsB);                                                                          // switch connectivities in lattice between 4 and 3 coordinate nodes
    void mixCnx34(VecF<int> mixIdsA, VecF<int> mixIdsB);                                                                                   // mix connectivities in lattice between 4 and 3 coordinate nodes
    void mixCnx(VecF<int> mixIdsA, VecF<int> mixIdsB);                                                                                     // mix connectivities in lattice between 4 and 3 coordinate nodes
    bool checkThreeRingEdges(int id);                                                                                                      // prevent edges being part of three rings
    bool convexRearrangement(int cnxType, VecF<int> switchIdsA, VecF<int> switchIdsB);                                                     // rearrange nodes after switchi to maintain convexity
    VecF<int> monteCarloSwitchMove(Network network, double &energy);                                                                       // monte carlo switching move
    VecF<int> SpiralmonteCarloSwitchMoveLAMMPS(int a, double &SimpleGrapheneEnergy, double &TersoffGrapheneEnergy, double &TriangleRaftEnergy, double &BilayerEnergy, double &BNEnergy, int Selected);
    VecF<int> monteCarloSwitchMoveLAMMPS(double &SimpleGrapheneEnergy, double &TersoffGrapheneEnergy, double &TriangleRaftEnergy, double &BilayerEnergy, double &BNEnergy, int Selected);
    VecF<int> monteCarloCostSwitchMove(double &cost, double &energy, double pTarget, double rTarget); // monte carlo switching move with cost function
    VecF<int> monteCarloMixMove(double &energy);                                                      // monte carlo mixing move
    double costFunction(double &pTarget, double &rTarget);                                            // cost function based on ring statistics and assortative mixing
    double globalPotentialEnergy(bool useIntx, bool restrict, Network network);                       // calculate potential energy of entire system
    //
    //    double globalPotentialEnergyTR(bool useIntx, bool restrict, Network network); //calculate potential energy of entire system
    //
    VecF<int> globalGeometryOptimisation(bool useIntx, bool restrict, Network network); // geometry optimise entire system
    //
    //    void globalGeometryOptimisationTR(bool useIntx, bool restrict); //geometry optimise entire system
    //
    VecF<int> localGeometryOptimisation(int centreA, int centreB, int extent, bool useIntx, bool restrict);                                    // geometry optimise subsection of system
    void generateHarmonics(int id, VecR<int> &bonds, VecR<double> &bondParams, VecR<int> &angles, VecR<double> &angleParams, Network network); // generate harmonic interactions

    void generateHarmonicsOnly(int id, VecR<int> &bonds, VecR<double> &bondParams, Network network); // generate harmonic interactions

    void generateRingIntersections(int rId, VecR<int> &intersections);   // generate ring intersection interactions
    void generateConvexIntersections(int nId, VecR<int> &intersections); // generate convex ring intersection interactions
    void wrapCoordinates();                                              // wrap coordinates if periodic
    void syncCoordinates();                                              // update geometry optimised coordinates to networks
    //
    void syncCoordinatesTD(); // update geometry optimised coordinates to networks
    //
    VecF<double> getNodeDistribution(string lattice);                                                    // get proportion of nodes of each size
    VecF<VecF<int>> getEdgeDistribution(string lattice);                                                 // get unnormalised proportion of node connections
    VecF<double> getAboavWeaire(string lattice);                                                         // get aboav-weaire parameters
    double getAssortativity(string lattice);                                                             // get network assortativity
    double getAboavWeaireEstimate(string lattice);                                                       // get estimate of aw alpha parameter from assortativity
    VecF<double> getEntropy(string lattice);                                                             // get node and edge distribution entropy
    VecF<double> getOptimisationGeometry(Network network, VecF<double> &lenHist, VecF<double> &angHist); // get bond/angle mean and standard deviation
    //
    VecF<double> getOptimisationGeometryTD(VecF<double> &lenHist, VecF<double> &angHist); // get bond/angle mean and standard deviation
    //
    void getRingAreas(VecF<double> &areaSum, VecF<double> &areaSqSum); // get sum of areas and squared areas of each ring size
    double getMaxCluster(string lattice, int nodeCnd);                 // get cluster statistics for given node coordination
    VecF<int> getMaxClusters(string lattice, int minCnd, int maxCnd);  // get cluster statistics for node coordinations
    bool checkConsistency();                                           // check networks are consistent
    bool checkCnxConsistency();                                        // check for mutual connections
    bool checkDescriptorConsistency();                                 // check descriptors are accurate
    bool checkConvexity();                                             // check all angles are convex
    bool checkConvexity(int id);                                       // check angles are convex around given node

    // Write Functions
    void writeXYZ(string prefix);
    void write(string prefix);
};

#endif // NL_LINKED_NETWORK_H
