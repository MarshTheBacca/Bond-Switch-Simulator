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
    // Member variables
    std::vector<Node> nodes;
    std::vector<double> dimensions;
    std::vector<double> reciprocalDimensions;
    std::string geometryCode;

    // Statistics
    double pearsonsCoeff;
    double entropy;
    std::map<int, double> nodeSizes;
    std::map<int, std::map<int, double>> assortativityDistribution;

    // Constructors
    Network();
    explicit Network(const int &numNodes);
    Network(const std::string &prefix, const LoggerPtr &logger); // construct by loading from files

    // Member Functions
    Network constructDual(const int &maxCnxs); // make dual graph
    void rescale(const double &scaleFactor);   // rescale coordinates

    void refreshStatistics();
    void refreshAssortativityDistribution();
    void refreshCoordinationDistribution();
    void refreshPearsonsCoeff();
    void refreshEntropy();

    double getAverageCoordination() const;
    double getAboavWeaire() const;
    double getAverageCoordination(const int &power) const;

    // Write functions
    void writeAux(std::ofstream &auxFile) const;
    void writeCrds(std::ofstream &crdFile) const;
    void writeCnxs(std::ofstream &cnxFile, const std::vector<std::vector<int>> &cnxs) const;
    std::vector<std::vector<int>> getNetCnxs() const;
    std::vector<std::vector<int>> getDualCnxs() const;
    void write(const std::string &prefix) const;

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
};

#endif // NL_NETWORK_H
