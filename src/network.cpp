#include "network.h"
#include "output_file.h"
#include "vector_tools.h"
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <ranges>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

const std::string BSS_NETWORK_PATH =
    std::filesystem::path("./input_files") / "bss_network";

/**
 * @brief Calculate the rounded square root of a number
 * @param num Number to calculate the rounded square root of
 * @return Rounded square root of the number
 */
inline int roundedSqrt(const int &num) {
  return static_cast<int>(std::round(std::sqrt(num)));
}

/**
 * @brief Convert a NetworkType to the corresponding networkString string
 * @param networkType NetworkType to convert
 * @return networkString string
 * @throw std::invalid_argument if networkType is not valid
 */
std::string NetworkTypeToString(NetworkType networkType) {
  switch (networkType) {
  case NetworkType::BASE_NETWORK:
    return "base_network";
  case NetworkType::DUAL_NETWORK:
    return "dual_network";
  default:
    throw std::invalid_argument("Invalid network type");
  }
}

/**
 * @brief Read a value from a file after a colon, used for reading the info file
 * @param infoFile File to read from
 * @return Value after the colon
 */
double readValueAfterColon(std::ifstream &infoFile) {
  std::string line;
  std::getline(infoFile, line);
  std::istringstream iss(line);
  std::string label;
  std::getline(iss, label, ':'); // Read until the colon
  double value;
  iss >> value; // Read the number after the colon
  return value;
}

/**
 * @brief default constructor
 */
Network::Network() = default;

/**
 * @brief Construct a network from files
 * @param networkString networkString of files to load
 * @param maxBaseCoordinationArg Maximum base coordination of nodes
 * @param maxDualCoordinationArg Maximum dual coordination of nodes
 * @param logger Logger to log to
 * @throw std::runtime_error if cannot open info file, crds file, net file or
 * dual file
 */
Network::Network(const NetworkType networkType, const LoggerPtr &logger)
    : type(networkType), networkString(NetworkTypeToString(networkType)) {
  logger->debug("Reading file: " + networkString + "_info.txt");
  readInfo(std::filesystem::path(BSS_NETWORK_PATH) /
           (networkString + "_info.txt"));
  nodes.reserve(numNodes);
  for (int i = 0; i < numNodes; ++i) {
    nodes.emplace_back(Node(i));
  }
  logger->debug("Reading file: " + networkString + "_coords.txt");
  readCoords(std::filesystem::path(BSS_NETWORK_PATH) /
             (networkString + "_coords.txt"));
  logger->debug("Reading file: " + networkString + "_connections.txt");
  readConnections(std::filesystem::path(BSS_NETWORK_PATH) /
                      (networkString + "_connections.txt"),
                  false);
  logger->debug("Reading file: " + networkString + "_dual_connections.txt");
  readConnections(std::filesystem::path(BSS_NETWORK_PATH) /
                      (networkString + "_dual_connections.txt"),
                  true);
}

/**
 * @Brief Read the number of nodes and dimensions from the info file
 * @param filePath Path to the info file
 * @param numNodes Number of nodes
 * @param dimensions Dimensions of the network
 * @throw std::runtime_error if info file not found
 */
void Network::readInfo(const std::string &filePath) {
  std::ifstream infoFile(filePath, std::ios::in);
  if (!infoFile.is_open()) {
    throw std::runtime_error("Cannot open info file: " + filePath);
  }
  // Read number of nodes
  numNodes = static_cast<int>(readValueAfterColon(infoFile));
  // Read dimensions
  dimensions[0] = readValueAfterColon(infoFile);
  dimensions[1] = readValueAfterColon(infoFile);
}

void Network::readCoords(const std::string &filePath) {
  std::ifstream coordsFile(filePath, std::ios::in);
  if (!coordsFile.is_open()) {
    throw std::runtime_error("Cannot coords open file: " + filePath);
  }
  std::array<double, 2> coord{0.0, 0.0};
  std::string line;
  std::ranges::for_each(nodes, [&coordsFile, &coord, &line](Node &node) {
    std::getline(coordsFile, line);
    std::istringstream ss(line);
    ss >> coord[0];
    ss >> coord[1];
    node.coord = coord;
  });
}

void Network::readConnections(const std::string &filePath, const bool &isDual) {
  std::ifstream connectionsFile(filePath, std::ios::in);
  if (!connectionsFile.is_open()) {
    throw std::runtime_error("Cannot open connections file: " + filePath);
  }
  int connectedNodeID;
  std::string line;
  std::ranges::for_each(
      nodes, [&connectionsFile, &connectedNodeID, &line, isDual](Node &node) {
        std::getline(connectionsFile, line);
        std::istringstream ss(line);
        while (ss >> connectedNodeID) {
          if (!isDual)
            node.netConnections.insert(connectedNodeID);
          else
            node.dualConnections.insert(connectedNodeID);
        }
      });
}

/**
 * @brief Find the number of unique dual nodes in all the node's dualConnections
 * @return Number of unique dual nodes
 */
int Network::findNumberOfUniqueDualNodes() {
  int numberOfUniqueDualNodes = -1;
  std::ranges::for_each(nodes, [&numberOfUniqueDualNodes](Node &node) {
    std::ranges::for_each(node.dualConnections,
                          [&numberOfUniqueDualNodes](int dualCnx) {
                            if (dualCnx > numberOfUniqueDualNodes)
                              numberOfUniqueDualNodes = dualCnx;
                          });
  });
  return numberOfUniqueDualNodes + 1;
}

/**
 * @brief Centres all nodes relative to their dual connections
 * @param pairedNetwork Base network to provide dual node IDs
 */
void Network::centreRings(const Network &pairedNetwork) {
  std::array<double, 2> coordSum{0.0, 0.0};
  // For every node
  std::ranges::for_each(nodes, [&pairedNetwork, &coordSum, this](Node &node) {
    // Reset coordSum to 0
    std::ranges::fill(coordSum, 0.0);

    // Get a sum of the relative vectors to each dual node
    std::ranges::for_each(node.dualConnections, [&pairedNetwork, &coordSum,
                                                 &node, this](int dualCnx) {
      std::array<double, 2> pbcCoords = pbcArray(
          node.coord, pairedNetwork.nodes[dualCnx].coord, this->dimensions);
      arrayAdd(coordSum, pbcCoords);
    });

    // Average the sum of the relative vectors
    divideArray(coordSum, static_cast<double>(node.numDualConnections()));

    // Move the node to this new average coordinate
    arrayAdd(node.coord, coordSum);

    // Wrap the new coordinates back into the box if they are outside
    std::ranges::for_each(node.coord, [this](double &coord) {
      while (coord < 0) {
        // While the coordinate is less than 0, keep adding the dimension
        coord += this->dimensions[0];
      }
      while (coord >= this->dimensions[0]) {
        // While the coordinate is ge the dimension, keep subtracting
        // the dimension
        coord -= this->dimensions[0];
      }
    });
  });
}

/**
 * @brief Rescale the network by a given factor
 * @param scaleFactor Factor to rescale by
 */
void Network::rescale(const double &scaleFactor) {
  containerMultiply(dimensions, scaleFactor);
  std::ranges::for_each(nodes, [&scaleFactor](Node &node) {
    containerMultiply(node.coord, scaleFactor);
  });
}

/**
 * @brief Refresh the entropy of node sizes
 */
void Network::refreshEntropy() {
  entropy = 0.0;
  std::ranges::for_each(nodeSizes, [this](const std::pair<int, double> &pair) {
    entropy -= pair.second * std::log(pair.second);
  });
}

/**
 * @brief Refresh the assortativity distribution of the network, which is a map
 * of maps
 */
void Network::refreshAssortativityDistribution() {
  assortativityDistribution.clear();
  std::ranges::for_each(nodes, [this](const Node &node) {
    std::ranges::for_each(node.netConnections, [&node, this](const int &cnx) {
      int cnxRingSize = nodes[cnx].numConnections();
      auto &outerMap = assortativityDistribution[node.numConnections()];
      if (!outerMap.contains(cnxRingSize)) {
        outerMap[cnxRingSize] = 1;
        return;
      }
      ++outerMap[cnxRingSize];
    });
  });
  std::ranges::for_each(assortativityDistribution,
                        [](std::pair<const int, std::map<int, double>> &pair) {
                          normaliseMap(pair.second);
                        });
}

/**
 * @brief Get the average coordination of the network
 * @return Average coordination
 */
double Network::getAverageCoordination() const {
  double totalCoordination = 0.0;
  std::ranges::for_each(nodes, [&totalCoordination](const Node &node) {
    totalCoordination += node.numConnections();
  });
  return totalCoordination / static_cast<double>(nodes.size());
}

/**
 * @brief Get the average coordination of the network to the power of a given
 * power, ie, <k^n>
 * @param power Power to raise the coordination to
 * @return Average coordination to the power of the given power
 */
double Network::getAverageCoordination(const int &power) const {
  double totalCoordination = 0.0;
  std::ranges::for_each(nodes, [&totalCoordination, power](const Node &node) {
    totalCoordination += std::pow(node.numConnections(), power);
  });
  return totalCoordination / static_cast<double>(nodes.size());
}

/**
 * @brief Refresh Pearson's correlation coefficient for the network

void Network::refreshPearsonsCoeff() {
    // r = <k>^2 / (<k><k^3> - <k^2>^2)  * sum_(j,k) ( jkP(j,k) - <k^2>^2)
    double avK = getAverageCoordination();                             // <k>
    double avKSquared = std::pow(avK, 2);                              // <k>^2
    double avKSquaredSquared = std::pow(getAverageCoordination(2), 2); //
<k^2>^2 double avKCubed = getAverageCoordination(3);                       //
<k^3> showNestedMap(assortativityDistribution); std::cout << "avK: " << avK << "
avKSquared: " << avKSquared << std::endl; std::cout << "avKSquaredSquared: " <<
avKSquaredSquared << " avKCubed: " << avKCubed << std::endl; double sum = 0.0;
    std::for_each(assortativityDistribution.begin(),
assortativityDistribution.end(), [&sum, avKSquaredSquared](const std::pair<int,
std::map<int, double>> &pair) { std::for_each(pair.second.begin(),
pair.second.end(), [&sum, avKSquaredSquared, &pair](const std::pair<int, double>
&innerPair) { sum += innerPair.first * pair.first * innerPair.second -
avKSquaredSquared;
        });
    });

    pearsonsCoeff = avKSquared / (avK * avKCubed - avKSquaredSquared) * sum;
    std::cout << "Pearsons Coeff: " << pearsonsCoeff << std::endl;
}
*/

void Network::refreshPearsonsCoeff() {
  double sum_x = 0;
  double sum_y = 0;
  double sum_xy = 0;
  double sum_x2 = 0;
  double sum_y2 = 0;
  int n = 0;

  for (const auto &outerPair : assortativityDistribution) {
    for (const auto &innerPair : outerPair.second) {
      double x = outerPair.first;
      double y = innerPair.first;
      double p = innerPair.second;

      sum_x += x * p;
      sum_y += y * p;
      sum_xy += x * y * p;
      sum_x2 += x * x * p;
      sum_y2 += y * y * p;
      n++;
    }
  }

  double mean_x = sum_x / n;
  double mean_y = sum_y / n;
  double cov_xy = sum_xy / n - mean_x * mean_y;
  double var_x = sum_x2 / n - mean_x * mean_x;
  double var_y = sum_y2 / n - mean_y * mean_y;

  pearsonsCoeff = cov_xy / (std::sqrt(var_x) * std::sqrt(var_y));
}

/**
 * @brief Refresh the probabilities of each node size
 */
void Network::refreshCoordinationDistribution() {
  nodeSizes.clear();
  std::ranges::for_each(nodes, [this](const Node &node) {
    try {
      ++nodeSizes.at(node.numConnections());
    } catch (std::out_of_range) {
      nodeSizes[node.numConnections()] = 1;
    }
  });
  normaliseMap(nodeSizes);
}

/**
 * @brief Write the info file with number of nodes and dimensions
 * @param infoFile File to write to
 */
void Network::writeInfo(std::ofstream &infoFile) const {
  infoFile << "Number of atoms: " << nodes.size() << "\n";
  if (!dimensions.empty()) {
    infoFile << "xhi: " << dimensions[0] << "\n";
    infoFile << "yhi: " << dimensions[1] << "\n";
  }
  infoFile.close();
}

/**
 * @brief Write coordinates to a file
 * @param crdFile File to write to
 */
void Network::writeCoords(std::ofstream &crdFile) const {
  crdFile << std::fixed << std::showpoint << std::setprecision(6);
  std::ranges::for_each(nodes, [&crdFile](const Node &node) {
    std::ranges::for_each(node.coord, [&crdFile](const double &coord) {
      crdFile << std::format("{:<20}", coord);
    });
    crdFile << std::endl;
  });
  crdFile.close();
}

/**
 * @brief Write connections to a file
 * @param connectionsFile File to write to
 * @param connections Connections to write
 */
void Network::writeConnections(
    std::ofstream &connectionsFile,
    const std::vector<std::set<int>> &connections) const {
  connectionsFile << std::fixed << std::showpoint << std::setprecision(1);
  for (int i = 0; i < nodes.size(); ++i) {
    std::ranges::for_each(
        connections[i], [&connectionsFile](const int &connection) {
          connectionsFile << std::format("{:<20}", connection);
        });
    connectionsFile << std::endl;
  }
  connectionsFile.close();
}

/**
 * @brief Get the net connections of the network
 * @return 2D vector of net connections
 */
std::vector<std::set<int>> Network::getConnections() const {
  std::vector<std::set<int>> netConnections(nodes.size());
  for (int i = 0; i < nodes.size(); ++i) {
    netConnections[i] = nodes[i].netConnections;
  }
  return netConnections;
}

/**
 * @brief Get the dual connections of the network
 * @return 2D vector of dual connections
 */
std::vector<std::set<int>> Network::getDualConnections() const {
  std::vector<std::set<int>> dualConnections{};
  std::ranges::for_each(nodes, [&dualConnections](const Node &node) {
    dualConnections.push_back(node.dualConnections);
  });
  return dualConnections;
}

/**
 * @brief Write network to files
 */
void Network::write() const {
  std::ofstream infoFile(std::filesystem::path("./output_files") /
                             (networkString + "_info.txt"),
                         std::ios::in | std::ios::trunc);
  writeInfo(infoFile);
  std::ofstream crdFile(std::filesystem::path("./output_files") /
                            (networkString + "_coords.txt"),
                        std::ios::in | std::ios::trunc);
  writeCoords(crdFile);
  std::ofstream netFile(std::filesystem::path("./output_files") /
                            (networkString + "_connections.txt"),
                        std::ios::in | std::ios::trunc);
  writeConnections(netFile, getConnections());
  std::ofstream dualFile(std::filesystem::path("./output_files") /
                             (networkString + "_dual_connections.txt"),
                         std::ios::in | std::ios::trunc);
  writeConnections(dualFile, getDualConnections());
}

/**
 * @brief Get the maximum coordination number for all nodes in the network
 * @return Maximum number of connections
 */
int Network::getMaxConnections() const {
  if (nodes.empty()) {
    throw std::runtime_error("Cannot get max connections of " + networkString +
                             ": no nodes");
  }
  auto maxNode =
      std::ranges::max_element(nodes, [](const Node &a, const Node &b) {
        return a.numConnections() < b.numConnections();
      });
  return maxNode->numConnections();
}

/**
 * @brief Get the maximum coordination number for all nodes in the network,
 * excluding those in excludeNodes
 * @param excludeNodes The IDs of all the nodes to exclude
 * @return Maximum number of connections
 */
int Network::getMaxConnections(const std::set<int> &excludeNodes) const {
  int maxConnections = 0;
  for (const auto &node : nodes) {
    if (excludeNodes.contains(node.id)) {
      continue;
    }
    if (node.numConnections() > maxConnections) {
      maxConnections = node.numConnections();
    }
  }
  return maxConnections;
}

/**
 * @brief Get the minimum number of connections of all nodes in the network
 * @return Minimum number of connections
 */
int Network::getMinConnections() const {
  auto minNode =
      std::ranges::min_element(nodes, [](const Node &a, const Node &b) {
        return a.netConnections.size() < b.netConnections.size();
      });
  return minNode->numConnections();
}

/**
 * @brief Get the minimum number of connections of all nodes in the network,
 * excluding those in excludeNodes
 * @param fixedNodes The IDs of all the fixed nodes
 * @return Minimum number of connections
 */
int Network::getMinConnections(const std::set<int> &excludeNodes) const {
  int minConnections = std::numeric_limits<int>::max();
  for (const auto &node : nodes) {
    if (excludeNodes.contains(node.id))
      continue;
    if (node.numConnections() < minConnections)
      minConnections = node.numConnections();
  }
  return minConnections;
}

/**
 * @brief Get the maximum dual coordination number for all nodes in the network
 * @return Maximum number of connections
 */
int Network::getMaxDualConnections() const {
  auto maxNode =
      std::ranges::max_element(nodes, [](const Node &a, const Node &b) {
        return a.numDualConnections() < b.numDualConnections();
      });
  return maxNode->numDualConnections();
}

/**
 * @brief Get the maximum dual coordination number for all nodes in the network,
 * excluding those in excludeNodes
 * @param excludeNodes The IDs of all the nodes to exclude
 * @return Maximum number of connections
 */
int Network::getMinDualConnections(const std::set<int> &excludeNodes) const {
  int minConnections = std::numeric_limits<int>::max();
  for (const auto &node : nodes) {
    if (excludeNodes.contains(node.id)) {
      continue;
    }
    if (node.numDualConnections() < minConnections) {
      minConnections = node.numDualConnections();
    }
  }
  return minConnections;
}

/**
 * @brief Get the minimum dual coordination number for all nodes in the network
 * @return Minimum number of connections
 */
int Network::getMinDualConnections() const {
  auto minNode =
      std::ranges::min_element(nodes, [](const Node &a, const Node &b) {
        return a.numDualConnections() < b.numDualConnections();
      });
  return minNode->numDualConnections();
}

/**
 * @brief Get the coordinates of the network in pairs
 * @return 1D vector of node coordinates
 */
std::vector<double> Network::getCoords() {
  std::vector<double> returnCoords;
  returnCoords.reserve(nodes.size() * 2);
  for (int i = 0; i < nodes.size(); i++) {
    returnCoords[i * 2] = nodes[i].coord[0];
    returnCoords[i * 2 + 1] = nodes[i].coord[1];
  }
  return returnCoords;
}

double Network::getAboavWeaire() const { return 0.0; }

void Network::refreshStatistics() {
  refreshCoordinationDistribution();
  refreshAssortativityDistribution();
  refreshPearsonsCoeff();
  refreshEntropy();
}

/**
 * @brief Display the network to the logger for debugging purposes
 * @param logger Logger to log to
 */
void Network::display(const LoggerPtr &logger) const {
  logger->info("Network type: {}", networkString);
  logger->info("Number of nodes: {}", numNodes);
  logger->info("Dimensions: [{}, {}]", dimensions[0], dimensions[1]);
  logger->info("Nodes:");
  std::ranges::for_each(
      nodes, [&logger](const Node &node) { logger->info(node.toString()); });
}