#include "network.h"
#include <filesystem>

const std::string BSS_NETWORK_PATH = std::filesystem::path("./input_files") / "bss_network";

/**
 * @brief Calculate the rounded square root of a number
 * @param num Number to calculate the rounded square root of
 * @return Rounded square root of the number
 */
inline int roundedSqrt(const int &num) {
    return std::round(std::sqrt(num));
}

/**
 * @brief Convert a NetworkType to the corresponding networkString string
 * @param networkType NetworkType to convert
 * @return networkString string
 * @throw std::invalid_argument if networkType is not valid
*/
std::string NetworkTypeToString(NetworkType networkType){
    switch(networkType){
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
double readValueAfterColon(std::ifstream& infoFile) {
    std::string line;
    std::getline(infoFile, line);
    std::istringstream iss(line);
    std::string label;
    std::getline(iss, label, ':');  // Read until the colon
    double value;
    iss >> value;  // Read the number after the colon
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
 * @throw std::runtime_error if cannot open info file, crds file, net file or dual file
 */
Network::Network(const NetworkType networkType, const LoggerPtr &logger) : type(networkType), networkString(NetworkTypeToString(networkType)) {
    logger->debug("Reading file: " + networkString + "_info.txt");
    readInfo(std::filesystem::path(BSS_NETWORK_PATH) / (networkString + "_info.txt"));
    nodes.reserve(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        Node node(i);
        nodes.push_back(node);
    }
    logger->debug("Reading file: " + networkString + "_coords.txt");
    readCoords(std::filesystem::path(BSS_NETWORK_PATH) / (networkString + "_coords.txt"));
    logger->debug("Reading file: " + networkString + "_connections.txt");
    readConnections(std::filesystem::path(BSS_NETWORK_PATH) / (networkString + "_connections.txt"), false);
    logger->debug("Reading file: " + networkString + "_dual_connections.txt");
    readConnections(std::filesystem::path(BSS_NETWORK_PATH) / (networkString + "_dual_connections.txt"), true);
}


/**
 * @Brief Read the number of nodes and dimensions from the info file
 * @param filePath Path to the info file
 * @param numNodes Number of nodes
 * @param dimensions Dimensions of the network
 * @throw std::runtime_error if info file not found
*/
void Network::readInfo(const std::string& filePath) {
    std::ifstream infoFile(filePath, std::ios::in);
    if (!infoFile.is_open()) {
        throw std::runtime_error("Cannot open info file: " + filePath);
    }
    // Read number of nodes
    numNodes = static_cast<int>(readValueAfterColon(infoFile));
    // Read dimensions
    dimensions = std::vector<double>(2);
    dimensions[0] = readValueAfterColon(infoFile);
    dimensions[1] = readValueAfterColon(infoFile);
}

void Network::readCoords(const std::string &filePath) {
    std::ifstream coordsFile(filePath, std::ios::in);
    if (!coordsFile.is_open()) {
        throw std::runtime_error("Cannot coords open file: " + filePath);
    }
    std::vector<double> crd(2);
    std::string line;
    std::for_each(nodes.begin(), nodes.end(), [&coordsFile, &crd, &line](Node &node) {
        std::getline(coordsFile, line);
        std::istringstream ss(line);
        ss >> crd[0];
        ss >> crd[1];
        node.crd = crd;
    });
}

void Network::readConnections(const std::string &filePath, const bool &isDual) {
    std::ifstream connectionsFile(filePath, std::ios::in);
    if (!connectionsFile.is_open()) {
        throw std::runtime_error("Cannot open connections file: " + filePath);
    }
    int cnx;
    std::string line;
    std::for_each(nodes.begin(), nodes.end(), [&connectionsFile, &cnx, &line, isDual](Node &node) {
        std::getline(connectionsFile, line);
        std::istringstream ss(line);
        while (ss >> cnx) {
            if (!isDual)
                node.netConnections.push_back(cnx);
            else
                node.dualConnections.push_back(cnx);
        }
    });
}

/**
 * @brief Find the number of unique dual nodes in all the node's dualConnections
 * @return Number of unique dual nodes
 */
int Network::findNumberOfUniqueDualNodes() {
    int numberOfUniqueDualNodes = -1;
    std::for_each(nodes.begin(), nodes.end(), [&numberOfUniqueDualNodes](Node &node) {
        std::for_each(node.dualConnections.begin(), node.dualConnections.end(), [&numberOfUniqueDualNodes](int dualCnx) {
            if (dualCnx > numberOfUniqueDualNodes)
                numberOfUniqueDualNodes = dualCnx;
        });
    });
    return numberOfUniqueDualNodes + 1;
}

/**
 * @brief Centres all nodes relative to their dual connections
 * @param baseNetwork Base network to provide dual node IDs
 */
void Network::centreRings(const Network &baseNetwork) {
    std::vector<double> total(2, 0.0);
    for (int ringNode = 0; ringNode < nodes.size(); ++ringNode) {
        Node &selectedRingNode = nodes[ringNode];
        total[0] = total[1] = 0.0;
        for (int neighbour = 0; neighbour < selectedRingNode.dualConnections.size(); ++neighbour) {
            std::vector<double> pbcCoords = pbcVector(selectedRingNode.crd,
                                                      baseNetwork.nodes[selectedRingNode.dualConnections[neighbour]].crd,
                                                      dimensions);
            total[0] += pbcCoords[0];
            total[1] += pbcCoords[1];
        }
        total[0] /= selectedRingNode.dualConnections.size();
        total[1] /= selectedRingNode.dualConnections.size();

        selectedRingNode.crd[0] += total[0];
        selectedRingNode.crd[1] += total[1];

        // Wrap the new coordinates back into the box
        for (size_t i = 0; i < selectedRingNode.crd.size(); ++i) {
            while (selectedRingNode.crd[i] < 0) {
                selectedRingNode.crd[i] += dimensions[i];
            }
            while (selectedRingNode.crd[i] >= dimensions[i]) {
                selectedRingNode.crd[i] -= dimensions[i];
            }
        }
    }
}

/**
 * @brief Rescale the network by a given factor
 * @param scaleFactor Factor to rescale by
 */
void Network::rescale(const double &scaleFactor) {
    vectorMultiply(dimensions, scaleFactor);
    std::for_each(nodes.begin(), nodes.end(), [&scaleFactor](Node &node) {
        vectorMultiply(node.crd, scaleFactor);
    });
}

/**
 * @brief Refresh the entropy of node sizes
 */
void Network::refreshEntropy() {
    entropy = 0.0;
    std::for_each(nodeSizes.begin(), nodeSizes.end(), [this](const std::pair<int, double> &pair) {
        entropy -= pair.second * std::log(pair.second);
    });
}

/**
 * @brief Refresh the assortativity distribution of the network, which is a map of maps
 */
void Network::refreshAssortativityDistribution() {
    assortativityDistribution.clear();
    std::for_each(nodes.begin(), nodes.end(), [this](const Node &node) {
        std::for_each(node.netConnections.begin(), node.netConnections.end(), [&node, this](const int &cnx) {
            int cnxRingSize = nodes[cnx].netConnections.size();
            auto &outerMap = assortativityDistribution[node.netConnections.size()];
            if (outerMap.find(cnxRingSize) == outerMap.end()) {
                outerMap[cnxRingSize] = 1;
            } else {
                ++outerMap[cnxRingSize];
            }
        });
    });
    std::for_each(assortativityDistribution.begin(), assortativityDistribution.end(), [](std::pair<const int, std::map<int, double>> &pair) {
        normaliseMap(pair.second);
    });
}

/**
 * @brief Get the average coordination of the network
 * @return Average coordination
 */
double Network::getAverageCoordination() const {
    double totalCoordination = 0.0;
    std::for_each(nodes.begin(), nodes.end(), [&totalCoordination](const Node &node) {
        totalCoordination += node.netConnections.size();
    });
    return totalCoordination / nodes.size();
}

/**
 * @brief Get the average coordination of the network to the power of a given power, ie, <k^n>
 * @param power Power to raise the coordination to
 * @return Average coordination to the power of the given power
 */
double Network::getAverageCoordination(const int &power) const {
    double totalCoordination = 0.0;
    std::for_each(nodes.begin(), nodes.end(), [&totalCoordination, power](const Node &node) {
        totalCoordination += std::pow(node.netConnections.size(), power);
    });
    return totalCoordination / nodes.size();
}

/**
 * @brief Refresh Pearson's correlation coefficient for the network

void Network::refreshPearsonsCoeff() {
    // r = <k>^2 / (<k><k^3> - <k^2>^2)  * sum_(j,k) ( jkP(j,k) - <k^2>^2)
    double avK = getAverageCoordination();                             // <k>
    double avKSquared = std::pow(avK, 2);                              // <k>^2
    double avKSquaredSquared = std::pow(getAverageCoordination(2), 2); // <k^2>^2
    double avKCubed = getAverageCoordination(3);                       // <k^3>
    showNestedMap(assortativityDistribution);
    std::cout << "avK: " << avK << " avKSquared: " << avKSquared << std::endl;
    std::cout << "avKSquaredSquared: " << avKSquaredSquared << " avKCubed: " << avKCubed << std::endl;
    double sum = 0.0;
    std::for_each(assortativityDistribution.begin(), assortativityDistribution.end(), [&sum, avKSquaredSquared](const std::pair<int, std::map<int, double>> &pair) {
        std::for_each(pair.second.begin(), pair.second.end(), [&sum, avKSquaredSquared, &pair](const std::pair<int, double> &innerPair) {
            sum += innerPair.first * pair.first * innerPair.second - avKSquaredSquared;
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
    std::for_each(nodes.begin(), nodes.end(), [this](const Node &node) {
        try {
            ++nodeSizes.at(node.netConnections.size());
        } catch (std::out_of_range &e) {
            nodeSizes[node.netConnections.size()] = 1;
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
        infoFile << "xhi: " << *dimensions.begin() << "\n";
        infoFile << "yhi: " << *(--dimensions.end()) << "\n";
    }
    infoFile.close();
}

/**
 * @brief Write coordinates to a file
 * @param crdFile File to write to
 */
void Network::writeCoords(std::ofstream &crdFile) const {
    crdFile << std::fixed << std::showpoint << std::setprecision(6);
    std::for_each(nodes.begin(), nodes.end(), [&crdFile](const Node &node) {
        std::for_each(node.crd.begin(), node.crd.end(), [&crdFile](const double &crd) {
            crdFile << std::setw(20) << std::left << crd;
        });
        crdFile << std::endl;
    });
    crdFile.close();
}

/**
 * @brief Write connections to a file
 * @param cnxFile File to write to
 * @param cnxs Connections to write
 */
void Network::writeConnections(std::ofstream &cnxFile, const std::vector<std::vector<int>> &cnxs) const {
    cnxFile << std::fixed << std::showpoint << std::setprecision(1);
    for (int i = 0; i < nodes.size(); ++i) {
        std::for_each(cnxs[i].begin(), cnxs[i].end(), [&cnxFile](const int &cnx) {
            cnxFile << std::setw(20) << std::left << cnx;
        });
        cnxFile << std::endl;
    }
    cnxFile.close();
}

/**
 * @brief Get the net connections of the network
 * @return 2D vector of net connections
 */
std::vector<std::vector<int>> Network::getConnections() const {
    std::vector<std::vector<int>> netConnections(nodes.size());
    for (int i = 0; i < nodes.size(); ++i) {
        netConnections[i] = nodes[i].netConnections;
    }
    return netConnections;
}

/**
 * @brief Get the dual connections of the network
 * @return 2D vector of dual connections
 */
std::vector<std::vector<int>> Network::getDualConnections() const {
    std::vector<std::vector<int>> dualConnections(nodes.size());
    for (int i = 0; i < nodes.size(); ++i) {
        dualConnections[i] = nodes[i].dualConnections;
    }
    return dualConnections;
}

/**
 * @brief Write network to files
 */
void Network::write() const {
    std::ofstream infoFile(std::filesystem::path("./output_files") / (networkString + "_info.txt"), std::ios::in | std::ios::trunc);
    writeInfo(infoFile);
    std::ofstream crdFile(std::filesystem::path("./output_files") / (networkString + "_coords.txt"), std::ios::in | std::ios::trunc);
    writeCoords(crdFile);
    std::ofstream netFile(std::filesystem::path("./output_files") / (networkString + "_connections.txt"), std::ios::in | std::ios::trunc);
    writeConnections(netFile, getConnections());
    std::ofstream dualFile(std::filesystem::path("./output_files") / (networkString + "_dual_connections.txt"), std::ios::in | std::ios::trunc);
    writeConnections(dualFile, getDualConnections());
}
/**
 * @brief Get the maximum coordination number for all nodes in the network
 * @return Maximum number of connections
 */
int Network::getMaxConnections() const {
    if (nodes.empty()) {
        throw std::runtime_error("Cannot get max connections of " + networkString + ": no nodes");
    }
    auto maxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.netConnections.size() < b.netConnections.size();
    });
    return maxNode->netConnections.size();
}

/**
 * @brief Get the maximum coordination number for all nodes in the network, excluding those in excludeNodes
 * @param excludeNodes The IDs of all the nodes to exclude
 * @return Maximum number of connections
 */
int Network::getMaxConnections(const std::unordered_set<int> &excludeNodes) const {
    int maxConnections = 0;
    for (const auto &node : nodes) {
        if (excludeNodes.find(node.id) != excludeNodes.end())
            continue;
        if (node.netConnections.size() > maxConnections)
            maxConnections = node.netConnections.size();
    }
    return maxConnections;
}

/**
 * @brief Get the minimum number of connections of all nodes in the network
 * @return Minimum number of connections
 */
int Network::getMinConnections() const {
    auto minNode = std::min_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.netConnections.size() < b.netConnections.size();
    });
    return minNode->netConnections.size();
}

/**
 * @brief Get the minimum number of connections of all nodes in the network, excluding those in excludeNodes
 * @param fixedNodes The IDs of all the fixed nodes
 * @return Minimum number of connections
 */
int Network::getMinConnections(const std::unordered_set<int> &excludeNodes) const {
    int minConnections = std::numeric_limits<int>::max();
    for (const auto &node : nodes) {
        if (excludeNodes.find(node.id) != excludeNodes.end())
            continue;
        if (node.netConnections.size() < minConnections)
            minConnections = node.netConnections.size();
    }
    return minConnections;
}

/**
 * @brief Get the maximum dual coordination number for all nodes in the network
 * @return Maximum number of connections
 */
int Network::getMaxDualConnections() const {
    auto maxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.dualConnections.size() < b.dualConnections.size();
    });
    return maxNode->dualConnections.size();
}

/**
 * @brief Get the maximum dual coordination number for all nodes in the network, excluding those in excludeNodes
 * @param excludeNodes The IDs of all the nodes to exclude
 * @return Maximum number of connections
 */
int Network::getMinDualConnections(const std::unordered_set<int> &excludeNodes) const {
    int minConnections = std::numeric_limits<int>::max();
    for (const auto &node : nodes) {
        if (excludeNodes.find(node.id) != excludeNodes.end())
            continue;
        if (node.dualConnections.size() < minConnections)
            minConnections = node.dualConnections.size();
    }
    return minConnections;
}

/**
 * @brief Get the minimum dual coordination number for all nodes in the network
 * @return Minimum number of connections
 */
int Network::getMinDualConnections() const {
    auto minNode = std::min_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.dualConnections.size() < b.dualConnections.size();
    });
    return minNode->dualConnections.size();
}

/**
 * @brief Get the coordinates of the network in pairs
 * @return 1D vector of node coordinates
 */
std::vector<double> Network::getCoords() {
    std::vector<double> returnCoords;
    returnCoords.reserve(nodes.size() * 2);
    for (int i = 0; i < nodes.size(); i++) {
        returnCoords[i * 2] = nodes[i].crd[0];
        returnCoords[i * 2 + 1] = nodes[i].crd[1];
    }
    return returnCoords;
}

double Network::getAboavWeaire() const {
    return 0.0;
}

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
void Network::display(const LoggerPtr &logger) const{
    logger->info("Network type: {}", networkString);
    logger->info("Number of nodes: {}", numNodes);
    logger->info("Dimensions: [{}, {}]", dimensions[0], dimensions[1]);
    logger->info("Nodes:");
    std::for_each(nodes.begin(), nodes.end(), [&logger](const Node &node) {
        logger->info(node.toString());
    });
}