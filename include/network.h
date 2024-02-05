// Network contains nodes and topological information
#ifndef NL_NETWORK_H
#define NL_NETWORK_H

#include "node.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>

#include <spdlog/spdlog.h>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

class Network {
public:
  // Data members
  VecF<double> pb;
  VecF<double> rpb; //(reciprocal) periodic boundary
  int maxNetCnxs;
  int maxDualCnxs;
  VecR<Node> nodes;                 // list of nodes
  std::string geometryCode;         // geometry of system
  VecF<int> nodeDistribution;       // number of each type of node
  VecF<VecF<int>> edgeDistribution; // number of each type of edge
  double SiScaling = 1.609 * sqrt(32.0 / 9.0) / 0.52917721090380;

  // Constructors
  Network();
  // Construct a triangular lattice
  explicit Network(int nNodes);
  Network(int nNodes, int maxCnxs);
  Network(const std::string &prefix, int maxNetCnxsA, int maxDualCnxsA,
          LoggerPtr logger); // construct by loading from files
  Network(VecR<Node> nodesA, VecF<double> pbA, VecF<double> rpbA,
          const std::string &type, int maxNetCnxsA, int maxDualCnxsA);

  // Member Functions
  Network constructDual(int maxCnxs, LoggerPtr logger); // make dual graph
  void gen2ndOrderConnections(Network dualNetwork);
  void rescale(double scaleFactor); // rescale coordinates
  void findLocalRegion(int a, int b, int extent, VecR<int> &local,
                       VecR<int> &fixedInner, VecR<int> &fixedOuter);
  VecF<double> getNodeDistribution();       // proportion of each type of node
  VecF<VecF<double>> getEdgeDistribution(); // proportion of each type of node
  VecF<double> aboavWeaireParams();         // calculate Aboav-Weaire parameters
  double assortativity();                   // calculate network assortativity
  double aboavWeaireEstimate();             // estimate aw alpha parameter
  VecF<double> entropy(); // calculate entropy of node and edge distribution
  VecF<int> maxClusters(
      int minCnd, int maxCnd, int minInnerCnxs,
      int minOuterCnxs); // get cluster statistics for node coordinations
  double
  maxCluster(int nodeCnd); // get cluster statistics for given node coordination
  void write(const std::string &prefix);
  void writeXYZ(const std::string &prefix, const std::string &element);
  void writeBilayerA(const std::string &prefix, float cutoff);
  void writeBilayerB(const std::string &prefix);
  void writeBN(const std::string &prefix);
  bool r_ij(int i, int j, float cutoff);
  int getMaxCnxs();
  int getMinCnxs();

private:
  // Refactored the initialiseTriangleLattice function
  void initialiseTriangularLattice(int dim);
  void assignCoordinates(int dim);
  void makeConnections(int dim, int dimSq);
  void addConnectionsBasedOnParity(int y, int id, int cnx, int dim, int dimSq);
  void addMoreConnectionsBasedOnParity(int y, int id, int cnx, int dim,
                                       int dimSq);
  void makeDualConnections(int dim, int dimSq);
  void addDualConnections(int y, int id, int dim, int dimSq2);
  int findNumberOfUniqueDualNodes();
  void addUnorderedDualConnections(Network &dualNetwork);
  void orderDualConnections(Network &dualNetwork);
  void addOrderedNetworkConnections(Network &dualNetwork);
  void setCoordinatesAtCentreOfDualConnections(Network &dualNetwork);

  void initialiseDescriptors(int maxCnxs);
};

#endif // NL_NETWORK_H
