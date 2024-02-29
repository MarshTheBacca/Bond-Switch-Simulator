#include "network.h"

/**
 * @brief Calculate the rounded square root of a number
 * @param num Number to calculate the rounded square root of
 * @return Rounded square root of the number
 */
inline int roundedSqrt(const int &num) {
    return std::round(std::sqrt(num));
}

/**
 * @brief default constructor
 */
Network::Network() = default;

/**
 * @brief Construct a triangular lattice
 * @param nNodes Number of nodes
 */
Network::Network(const int &nNodes) {
    initialiseTriangularLattice(roundedSqrt(nNodes));
}

/**
 * @brief Construct a network from files
 * @param prefix Prefix of files to load
 * @param maxBaseCoordinationArg Maximum base coordination of nodes
 * @param maxDualCoordinationArg Maximum dual coordination of nodes
 * @param logger Logger to log to
 * @throw std::runtime_error if cannot open aux file, crds file, net file or dual file
 */
Network::Network(const std::string &pathPrefix, const LoggerPtr &logger) {
    logger->debug("Reading aux file {} ...", pathPrefix + "_aux.dat");

    // Initialise variables with aux file information
    std::istringstream ss("");
    std::ifstream auxFile(pathPrefix + "_aux.dat", std::ios::in);
    if (!auxFile.is_open()) {
        throw std::runtime_error("Aux file not found!");
    }

    int numNodes;
    std::string line;
    std::getline(auxFile, line);
    std::istringstream(line) >> numNodes;
    logger->debug("Number of nodes: {}", numNodes);
    std::getline(auxFile, line);

    std::getline(auxFile, line);
    std::istringstream(line) >> geometryCode;
    dimensions = std::vector<double>(2);
    reciprocalDimensions = std::vector<double>(2);
    std::getline(auxFile, line);
    ss.str(line);
    ss >> dimensions[0];
    ss >> dimensions[1];
    std::getline(auxFile, line);
    ss.str(line);
    ss >> reciprocalDimensions[0];
    ss >> reciprocalDimensions[1];
    auxFile.close();

    nodes.reserve(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        Node node(i);
        nodes.push_back(node);
    }

    logger->debug("Reading crds file {} ...", pathPrefix + "_crds.dat");
    std::ifstream crdFile(pathPrefix + "_crds.dat", std::ios::in);
    if (!crdFile.is_open()) {
        throw std::runtime_error("crds.dat file not found!");
    }
    std::vector<double> crd(2);

    std::for_each(nodes.begin(), nodes.end(), [&crdFile, &crd, &line](Node &node) {
        std::getline(crdFile, line);
        std::istringstream ss(line);
        ss >> crd[0];
        ss >> crd[1];
        node.crd = crd;
    });

    crdFile.close();

    logger->debug("Reading net file {} ...", pathPrefix + "_net.dat");
    std::ifstream netFile(pathPrefix + "_net.dat", std::ios::in);
    if (!netFile.is_open()) {
        throw std::runtime_error("net.dat file not found!");
    }
    int cnx;
    std::for_each(nodes.begin(), nodes.end(), [&netFile, &cnx, &line](Node &node) {
        std::getline(netFile, line);
        std::istringstream ss(line);
        while (ss >> cnx) {
            node.netCnxs.push_back(cnx);
        }
    });

    netFile.close();

    // Read dual connections
    logger->debug("Reading dual file {} ...", pathPrefix + "_dual.dat");
    std::ifstream dualFile(pathPrefix + "_dual.dat", std::ios::in);
    if (!dualFile.is_open()) {
        throw std::runtime_error("dual.dat file not found!");
    }
    std::for_each(nodes.begin(), nodes.end(), [&dualFile, &cnx, &line](Node &node) {
        std::getline(dualFile, line);
        std::istringstream ss(line);
        while (ss >> cnx) {
            node.dualCnxs.push_back(cnx);
        }
    });

    dualFile.close();

    // Set up descriptors
    logger->debug("Number of nodes: {}", nodes.size());
}

// Initialise triangular lattice of periodic 6-coordinate nodes
void Network::initialiseTriangularLattice(const int &dim) {
    geometryCode = "2DE"; // 2D euclidean
    int dimSq = dim * dim;

    // make 6 coordinate nodes
    nodes.reserve(dimSq);
    for (int i = 0; i < nodes.size(); ++i) {
        Node node(i);
        nodes.push_back(node);
    }

    // assign coordinates in layers, with unit bond lengths
    dimensions.reserve(2);
    reciprocalDimensions.reserve(2);
    dimensions[0] = dim;
    dimensions[1] = dim * sqrt(3) * 0.5;
    reciprocalDimensions[0] = 1.0 / dimensions[0];
    reciprocalDimensions[1] = 1.0 / dimensions[1];
    assignCoordinates(dim);
    makeConnections(dim, dimSq);
    makeDualConnections(dim, dimSq);
}

void Network::assignCoordinates(const int &dim) {
    std::vector<double> c(2);
    double dy = sqrt(3.0) * 0.5;
    for (int y = 0; y < dim; ++y) {
        c[1] = 0.5 * dy + y * dy;
        for (int x = 0; x < dim; ++x) {
            c[0] = 0.5 * (y % 2) + x;
            nodes[x + y * dim].crd = c;
        }
    }
}

void Network::makeConnections(const int &dim, const int &dimSq) {
    int id = 0;
    int cnx;
    for (int y = 0; y < dim; ++y) {
        for (int x = 0; x < dim; ++x) {
            cnx = y * dim + (id + dim - 1) % dim;
            nodes[id].netCnxs.push_back(cnx);
            addConnectionsBasedOnParity(y, id, cnx, dim, dimSq);
            addMoreConnectionsBasedOnParity(y, id, cnx, dim, dimSq);
            ++id;
        }
    }
}

void Network::addConnectionsBasedOnParity(const int &y, const int &id, int &cnx, const int &dim,
                                          const int &dimSq) {
    if (y % 2 == 0) {
        cnx = (cnx + dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
        cnx = (id + dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
    } else {
        cnx = (id + dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
        cnx = ((y + 1) * dim + (id + 1) % dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
    }
    cnx = y * dim + (id + 1) % dim;
    nodes[id].netCnxs.push_back(cnx);
}

void Network::addMoreConnectionsBasedOnParity(const int &y, const int &id, int &cnx, const int &dim,
                                              const int &dimSq) {
    if (y % 2 == 0) {
        cnx = (id + dimSq - dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
        cnx = (dimSq + (y - 1) * dim + (id + dim - 1) % dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
    } else {
        cnx = (dimSq + (y - 1) * dim + (id + dim + 1) % dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
        cnx = (dimSq + (y - 1) * dim + (id + dim) % dim) % dimSq;
        nodes[id].netCnxs.push_back(cnx);
    }
}

void Network::makeDualConnections(const int &dim, const int &dimSq) {
    int id = 0;
    int dimSq2 = 2 * dimSq;
    for (int y = 0; y < dim; ++y) {
        for (int x = 0; x < dim; ++x) {
            addDualConnections(y, id, dim, dimSq2);
            ++id;
        }
    }
}

void Network::addDualConnections(const int &y, const int &id, const int &dim, const int &dimSq2) {
    int cnx = (2 * y * dim + (id + dim - 1) % dim);
    nodes[id].dualCnxs.push_back(cnx);
    cnx = ((2 * y + 1) * dim + id % dim);
    nodes[id].dualCnxs.push_back(cnx);
    cnx = (2 * y * dim + id % dim);
    nodes[id].dualCnxs.push_back(cnx);

    if (y % 2 == 0) {
        cnx = (2 * y * dim + id % dim - dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.push_back(cnx);
        cnx = (2 * (y - 1) * dim + (id + dim - 1) % dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.push_back(cnx);
        cnx = (2 * (y - 1) * dim + (id + dim - 1) % dim + dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.push_back(cnx);
    } else {
        cnx = (2 * y * dim + (id + 1) % dim - dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.push_back(cnx);
        cnx = (2 * (y - 1) * dim + id % dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.push_back(cnx);
        cnx = (2 * (y - 1) * dim + id % dim + dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.push_back(cnx);
    }
}

Network Network::constructDual(const int &maxCnxs) {
    int numberOfUniqueDualNodes = findNumberOfUniqueDualNodes();
    Network dualNetwork(numberOfUniqueDualNodes);

    addUnorderedDualConnections(dualNetwork);
    orderDualConnections(dualNetwork);
    addOrderedNetworkConnections(dualNetwork);
    centreRings(dualNetwork);

    // set remaining parameters
    dualNetwork.dimensions = dimensions;
    dualNetwork.reciprocalDimensions = reciprocalDimensions;
    dualNetwork.geometryCode = geometryCode;

    return dualNetwork;
}
/**
 * @brief Find the number of unique dual nodes in all the node's dualCnxs
 * @return Number of unique dual nodes
 */
int Network::findNumberOfUniqueDualNodes() {
    int numberOfUniqueDualNodes = -1;
    std::for_each(nodes.begin(), nodes.end(), [&numberOfUniqueDualNodes](Node &node) {
        std::for_each(node.dualCnxs.begin(), node.dualCnxs.end(), [&numberOfUniqueDualNodes](int dualCnx) {
            if (dualCnx > numberOfUniqueDualNodes)
                numberOfUniqueDualNodes = dualCnx;
        });
    });
    return numberOfUniqueDualNodes + 1;
}

void Network::addUnorderedDualConnections(Network &dualNetwork) {
    for (int i = 0; i < nodes.size(); ++i) {
        std::for_each(nodes[i].dualCnxs.begin(), nodes[i].dualCnxs.end(),
                      [&dualNetwork, i](int dualConnectionId) {
                          dualNetwork.nodes[dualConnectionId].dualCnxs.push_back(i);
                      });
    }
}

void Network::orderDualConnections(Network &dualNetwork) {
    for (int i = 0; i < dualNetwork.nodes.size(); ++i) {
        std::vector<int> initialDualCnxs = dualNetwork.nodes[i].dualCnxs;
        std::vector<int> orderedCnxs;
        orderedCnxs.push_back(initialDualCnxs[0]);
        for (int j = 1; j < initialDualCnxs.size(); ++j) {
            std::unordered_set<int> common = intersectVectors(nodes[orderedCnxs[j - 1]].netCnxs, initialDualCnxs);
            auto it = common.begin();
            if (!vectorContains(orderedCnxs, *it))
                orderedCnxs.push_back(*it);
            else {
                ++it;
                orderedCnxs.push_back(*it);
            }
        }
        dualNetwork.nodes[i].dualCnxs = orderedCnxs;
    }
}

void Network::addOrderedNetworkConnections(Network &dualNetwork) {
    for (int i = 0; i < dualNetwork.nodes.size(); ++i) {
        std::vector<int> dualCnxs = dualNetwork.nodes[i].dualCnxs;
        std::unordered_set<int> common;
        for (int j = 0; j < dualCnxs.size(); ++j) {
            int k = (j + 1) % dualCnxs.size();
            common = intersectVectors(nodes[dualCnxs[j]].dualCnxs,
                                      nodes[dualCnxs[k]].dualCnxs);
            common.erase(i);
            dualNetwork.nodes[i].netCnxs.push_back(*common.begin());
        }
    }
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
        for (int neighbour = 0; neighbour < selectedRingNode.dualCnxs.size(); ++neighbour) {
            std::vector<double> pbcCoords = pbcVector(selectedRingNode.crd,
                                                      baseNetwork.nodes[selectedRingNode.dualCnxs[neighbour]].crd,
                                                      dimensions);
            total[0] += pbcCoords[0];
            total[1] += pbcCoords[1];
        }
        total[0] /= selectedRingNode.dualCnxs.size();
        total[1] /= selectedRingNode.dualCnxs.size();

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
    vectorDivide(reciprocalDimensions, scaleFactor);
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
        std::for_each(node.netCnxs.begin(), node.netCnxs.end(), [&node, this](const int &cnx) {
            int cnxRingSize = nodes[cnx].netCnxs.size();
            auto &outerMap = assortativityDistribution[node.netCnxs.size()];
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
        totalCoordination += node.netCnxs.size();
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
        totalCoordination += std::pow(node.netCnxs.size(), power);
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
            ++nodeSizes.at(node.netCnxs.size());
        } catch (std::out_of_range &e) {
            nodeSizes[node.netCnxs.size()] = 1;
        }
    });
    normaliseMap(nodeSizes);
}

/**
 * @brief Write the aux file with number of nodes, connectivity limits, geometry code and dimensions
 * @param auxFile File to write to
 */
void Network::writeAux(std::ofstream &auxFile) const {
    auxFile << std::fixed << std::showpoint << std::setprecision(1);
    auxFile << std::setw(10) << std::left << nodes.size() << std::endl;
    auxFile << std::setw(10) << std::left << "1"
            << std::setw(10) << std::left << "1" << std::endl;
    auxFile << std::setw(10) << std::left << geometryCode << std::endl;
    auxFile << std::fixed << std::showpoint << std::setprecision(6);
    std::for_each(dimensions.begin(), dimensions.end(), [&auxFile](const double &dimension) {
        auxFile << std::setw(20) << std::left << dimension;
    });
    auxFile << std::endl;
    std::for_each(reciprocalDimensions.begin(), reciprocalDimensions.end(), [&auxFile](const double &reciprocalDimension) {
        auxFile << std::setw(20) << std::left << reciprocalDimension;
    });
    auxFile << std::endl;
    auxFile.close();
}

/**
 * @brief Write coordinates to a file
 * @param crdFile File to write to
 */
void Network::writeCrds(std::ofstream &crdFile) const {
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
void Network::writeCnxs(std::ofstream &cnxFile, const std::vector<std::vector<int>> &cnxs) const {
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
std::vector<std::vector<int>> Network::getNetCnxs() const {
    std::vector<std::vector<int>> netCnxs(nodes.size());
    for (int i = 0; i < nodes.size(); ++i) {
        netCnxs[i] = nodes[i].netCnxs;
    }
    return netCnxs;
}

/**
 * @brief Get the dual connections of the network
 * @return 2D vector of dual connections
 */
std::vector<std::vector<int>> Network::getDualCnxs() const {
    std::vector<std::vector<int>> dualCnxs(nodes.size());
    for (int i = 0; i < nodes.size(); ++i) {
        dualCnxs[i] = nodes[i].dualCnxs;
    }
    return dualCnxs;
}

/**
 * @brief Write network to files
 * @param prefix The path and prefix of files to write
 */
void Network::write(const std::string &prefix) const {
    std::ofstream auxFile(prefix + "_aux.dat", std::ios::in | std::ios::trunc);
    writeAux(auxFile);
    std::ofstream crdFile(prefix + "_crds.dat", std::ios::in | std::ios::trunc);
    writeCrds(crdFile);
    std::ofstream netFile(prefix + "_net.dat", std::ios::in | std::ios::trunc);
    writeCnxs(netFile, getNetCnxs());
    std::ofstream dualFile(prefix + "_dual.dat", std::ios::in | std::ios::trunc);
    writeCnxs(dualFile, getDualCnxs());
}
/**
 * @brief Get the maximum coordination number for all nodes in the network
 * @return Maximum number of connections
 */
int Network::getMaxCnxs() const {
    auto maxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.netCnxs.size() < b.netCnxs.size();
    });
    return maxNode->netCnxs.size();
}

/**
 * @brief Get the maximum coordination number for all nodes in the network, excluding those in excludeNodes
 * @param excludeNodes The IDs of all the nodes to exclude
 * @return Maximum number of connections
 */
int Network::getMaxCnxs(const std::unordered_set<int> &excludeNodes) const {
    int maxCnxs = 0;
    for (const auto &node : nodes) {
        if (excludeNodes.find(node.id) != excludeNodes.end())
            continue;
        if (node.netCnxs.size() > maxCnxs)
            maxCnxs = node.netCnxs.size();
    }
    return maxCnxs;
}

/**
 * @brief Get the minimum number of connections of all nodes in the network
 * @return Minimum number of connections
 */
int Network::getMinCnxs() const {
    auto minNode = std::min_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.netCnxs.size() < b.netCnxs.size();
    });
    return minNode->netCnxs.size();
}

/**
 * @brief Get the minimum number of connections of all nodes in the network, excluding those in excludeNodes
 * @param fixedNodes The IDs of all the fixed nodes
 * @return Minimum number of connections
 */
int Network::getMinCnxs(const std::unordered_set<int> &excludeNodes) const {
    int minCnxs = std::numeric_limits<int>::max();
    for (const auto &node : nodes) {
        if (excludeNodes.find(node.id) != excludeNodes.end())
            continue;
        if (node.netCnxs.size() < minCnxs)
            minCnxs = node.netCnxs.size();
    }
    return minCnxs;
}

/**
 * @brief Get the maximum dual coordination number for all nodes in the network
 * @return Maximum number of connections
 */
int Network::getMaxDualCnxs() const {
    auto maxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.dualCnxs.size() < b.dualCnxs.size();
    });
    return maxNode->dualCnxs.size();
}

/**
 * @brief Get the maximum dual coordination number for all nodes in the network, excluding those in excludeNodes
 * @param excludeNodes The IDs of all the nodes to exclude
 * @return Maximum number of connections
 */
int Network::getMinDualCnxs(const std::unordered_set<int> &excludeNodes) const {
    int minCnxs = std::numeric_limits<int>::max();
    for (const auto &node : nodes) {
        if (excludeNodes.find(node.id) != excludeNodes.end())
            continue;
        if (node.dualCnxs.size() < minCnxs)
            minCnxs = node.dualCnxs.size();
    }
    return minCnxs;
}

/**
 * @brief Get the minimum dual coordination number for all nodes in the network
 * @return Minimum number of connections
 */
int Network::getMinDualCnxs() const {
    auto minNode = std::min_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.dualCnxs.size() < b.dualCnxs.size();
    });
    return minNode->dualCnxs.size();
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
