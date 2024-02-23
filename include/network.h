// Network contains nodes and topological information
#ifndef NL_NETWORK_H
#define NL_NETWORK_H

#include "node.h"
#include "vector_tools.h"
#include <cmath>
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
    int maxNetCnxs; // Maximum coordination number of nodes to nodes
    int maxDualCnxs;
    std::vector<double> dimensions;                 // Periodic boundary of network, xlo = ylo = 0, so dimensions = [xhi, yhi]
    std::vector<double> rpb;                        // Reciprical periodic boundary = [1/xhi, 1/yhi]
    std::string geometryCode;                       // geometry of system
                                                    // Maximum coordination number of nodes to dual nodes
    std::vector<Node> nodes;                        // list of nodes
    std::vector<int> nodeDistribution;              // number of each type of node
    std::vector<std::vector<int>> edgeDistribution; // number of each type of edge

    // Construct a triangular lattice
    explicit Network(const int &nNodes);
    Network();

    Network(const int &nNodes, const int &maxCnxs);
    Network(const std::string &prefix, const int &maxNetCnxsA, const int &maxDualCnxsA, LoggerPtr logger); // construct by loading from files

    // Member Functions
    Network constructDual(const int &maxCnxs); // make dual graph
    void rescale(const double &scaleFactor);   // rescale coordinates

    std::vector<double> getNodeDistribution() const;
    std::vector<std::vector<double>> getEdgeDistribution() const;    // proportion of each type of node
    std::tuple<double, double, double> getAboavWeaireParams() const; // calculate Aboav-Weaire parameters
    double getAssortativity();                                       // calculate network assortativity
    double getAboavWeaireEstimate();                                 // estimate aw alpha parameter
    std::vector<double> getEntropy();                                // calculate entropy of node and edge distribution

    void write(const std::string &prefix);
    bool r_ij(int i, int j, double cutoff);

    int getMaxCnxs();
    int getMaxDualCnxs();
    int getMinCnxs();

    std::vector<double> getCoords();
    void getCoords(std::vector<double> &coords);

  private:
    // Refactored the initialiseTriangleLattice function
    void initialiseTriangularLattice(const int &dim);
    void assignCoordinates(const int &dim);
    void makeConnections(const int &dim, const int &dimSq);
    void addConnectionsBasedOnParity(const int &y, const int &id, int &cnx, const int &dim, const int &dimSq);
    void addMoreConnectionsBasedOnParity(const int &y, const int &id, int &cnx, const int &dim, const int &dimSq);
    void makeDualConnections(const int &dim, const int &dimSq);
    void addDualConnections(const int &y, const int &id, const int &dim, const int &dimSq2);
    int findNumberOfUniqueDualNodes();
    void addUnorderedDualConnections(Network &dualNetwork);
    void orderDualConnections(Network &dualNetwork);
    void addOrderedNetworkConnections(Network &dualNetwork);
    void setCoordinatesAtCentreOfDualConnections(Network &dualNetwork);

    void initialiseDescriptors(const int &maxCnxs);
};

#endif // NL_NETWORK_H
