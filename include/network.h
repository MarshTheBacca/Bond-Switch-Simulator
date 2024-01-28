// Network contains nodes and topological information
#ifndef NL_NETWORK_H
#define NL_NETWORK_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>
#include <memory>
#include "node.h"

#include <spdlog/spdlog.h>

using namespace std;
using LoggerPtr = std::shared_ptr<spdlog::logger>;

class Network
{

private:
    // Default lattices
    void initialiseSquareLattice(int dim, int &maxCnxs);                        // 4 coordinate nodes, forming periodic square lattice
    void initialiseTriangularLattice(int dim, int &maxCnxs);                    // 6 coordinate nodes, forming periodic square lattice
    void initialiseAltSquareLattice(int dim, int &maxCnxs);                     // 4/2 coordinate nodes, forming periodic square lattice
    void initialiseMixedTSLattice(int dim, int &maxCnxs, double mixProportion); // 4&6 coordinate nodes, forming mixed triangular and square lattice
    void initialiseCubicLattice(int nNodes, int &maxCnxs);                      // 3&4 coordinate nodes, forming cubic lattice
    void initialiseGeodesicLattice(int nNodes, int &maxCnxs);                   // 5&6 coordinate nodes, forming geodesic lattice
    void initialiseDescriptors(int maxCnxs, LoggerPtr logger);                  // node descriptors

public:
    // Data members
    VecF<double> pb, rpb; //(reciprocal) periodic boundary
    int maxNetCnxs, maxDualCnxs;
    VecR<Node> nodes;                 // list of nodes
    string geometryCode;              // geometry of system
    VecF<int> nodeDistribution;       // number of each type of node
    VecF<VecF<int>> edgeDistribution; // number of each type of edge
    double SiScaling = 1.609 * sqrt((32.0 / 9.0)) / 0.52917721090380;

    // Constructors
    Network();
    Network(int nNodes, int maxCnxs);
    Network(int nNodes, string lattice, int maxCnxs, LoggerPtr logger);          // construct with default lattice
    Network(string prefix, int maxNetCnxsA, int maxDualCnxsA, LoggerPtr logger); // construct by loading from files
    Network(VecR<Node> nodesA, VecF<double> pbA, VecF<double> rpbA, string type, int maxNetCnxsA, int maxDualCnxsA);

    // Member Functions
    Network constructDual(int maxCnxs, LoggerPtr logger);          // make dual graph
    void generateAuxConnections(Network dualNetwork, int auxType); // generate auxilary connections
    void rescale(double scaleFactor);                              // rescale coordinates
    void project(string projType, double param);                   // project lattice onto different geometry
    void findLocalRegion(int a, int b, int extent, VecR<int> &local, VecR<int> &fixedInner, VecR<int> &fixedOuter);
    VecF<double> getNodeDistribution();                                                // proportion of each type of node
    VecF<VecF<double>> getEdgeDistribution();                                          // proportion of each type of node
    VecF<double> aboavWeaireParams();                                                  // calculate Aboav-Weaire parameters
    double assortativity();                                                            // calculate network assortativity
    double aboavWeaireEstimate();                                                      // estimate aw alpha parameter
    VecF<double> entropy();                                                            // calculate entropy of node and edge distribution
    VecF<int> maxClusters(int minCnd, int maxCnd, int minInnerCnxs, int minOuterCnxs); // get cluster statistics for node coordinations
    double maxCluster(int nodeCnd);                                                    // get cluster statistics for given node coordination
    void write(string prefix);
    void writeXYZ(string prefix, string element);
    void writeBilayerA(string prefix, float cutoff);
    void writeBilayerB(string prefix);
    void writeBN(string prefix);
    bool r_ij(int i, int j, float cutoff);
};

#endif // NL_NETWORK_H
