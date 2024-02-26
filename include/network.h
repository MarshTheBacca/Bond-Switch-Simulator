// Network contains nodes and topological information
#ifndef NL_NETWORK_H
#define NL_NETWORK_H

#include "node.h"
#include "output_file.h"
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
    int minNetCnxs; // Minimum coordination number of nodes to nodes
    int maxNetCnxs; // Maximum coordination number of nodes to nodes
    int minDualCnxs;
    int maxDualCnxs;
    std::vector<double> dimensions;                 // Periodic boundary of network, xlo = ylo = 0, so dimensions = [xhi, yhi]
    std::vector<double> reciprocalDimensions;       // Reciprical periodic boundary = [1/xhi, 1/yhi]
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
    double getAssortativityOld() const;
    double getAboavWeaireEstimate() const;
    std::vector<double> getEntropyOld() const; // calculate entropy of node and edge distribution

    double getAssortativity() const;
    double getAboavWeaire() const;
    double getEntropy() const;

    std::vector<int> getCoordinations(const int &minRingSize, const int &maxRingSize) const;

    // Write functions
    void writeAux(std::ofstream &auxFile) const;
    void writeCrds(std::ofstream &crdFile) const;
    void writeCnxs(std::ofstream &cnxFile, const std::vector<std::vector<int>> &cnxs) const;
    std::vector<std::vector<int>> getNetCnxs() const;
    std::vector<std::vector<int>> getDualCnxs() const;
    void write(const std::string &prefix) const;

    bool r_ij(int i, int j, double cutoff);

    int getMaxCnxs() const;
    int getMaxCnxs(const std::unordered_set<int> &fixedNodes) const;

    int getMinCnxs() const;
    int getMinCnxs(const std::unordered_set<int> &fixedNodes) const;

    int getMaxDualCnxs() const;
    int getMinDualCnxs() const;
    int getMinDualCnxs(const std::unordered_set<int> &fixedNodes) const;

    std::vector<double> getCoords();
    void centreRings(const Network &baseNetwork);

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

    void initialiseDescriptors(const int &maxCnxs);
};

#endif // NL_NETWORK_H
