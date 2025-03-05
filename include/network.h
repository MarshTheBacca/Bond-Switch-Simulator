// Network contains nodes and topological information
#ifndef NL_NETWORK_H
#define NL_NETWORK_H

#include "node.h"
#include "types.h"
#include <array>
#include <map>
#include <set>
#include <spdlog/spdlog.h>
#include <string>
#include <unordered_set>
#include <vector>

enum class NetworkType { BASE_NETWORK, DUAL_NETWORK };

struct Network {
  // Member variables
  NetworkType type;
  std::string networkString;
  int numNodes;
  std::vector<Node> nodes;
  std::array<double, 2> dimensions;

  // Statistics
  double pearsonsCoeff;
  double entropy;
  // Map of node degree to probability
  std::map<int, double> nodeSizes;
  std::map<int, std::map<int, double>> assortativityDistribution;

  // Constructors
  Network();
  Network(const NetworkType networkType,
          const LoggerPtr &logger); // construct by loading from files
  void readInfo(const std::string &filePath);
  void readCoords(const std::string &filePath);
  void readConnections(const std::string &filePath, const bool &isDual);

  // Member Functions
  void rescale(const double &scaleFactor); // rescale coordinates
  void refreshStatistics();
  void refreshAssortativityDistribution();
  void refreshCoordinationDistribution();
  void refreshPearsonsCoeff();
  void refreshEntropy();

  double getAverageCoordination() const;
  double getAboavWeaire() const;
  double getAverageCoordination(const int power) const;

  // Write functions
  void writeInfo(std::ofstream &infoFile) const;
  void writeCoords(std::ofstream &crdFile) const;
  void writeConnections(std::ofstream &cnxFile,
                        const std::vector<std::set<uint16_t>> &cnxs) const;
  std::vector<std::set<uint16_t>> getConnections() const;
  std::vector<std::set<uint16_t>> getDualConnections() const;
  void write() const;

  size_t getMaxConnections() const;
  size_t
  getMaxConnections(const std::unordered_set<uint16_t> &excludeNodes) const;

  size_t getMinConnections() const;
  size_t
  getMinConnections(const std::unordered_set<uint16_t> &excludeNodes) const;

  size_t getMaxDualConnections() const;
  size_t getMinDualConnections() const;
  size_t
  getMinDualConnections(const std::unordered_set<uint16_t> &excludeNodes) const;

  std::vector<double> getCoords();
  void centreRings(const Network &baseNetwork);

  int findNumberOfUniqueDualNodes();
  void display(const LoggerPtr &logger) const;
  Node &getRandomNode();
  Node &getRandomNodeConnection(const Node &node);
};

#endif // NL_NETWORK_H
