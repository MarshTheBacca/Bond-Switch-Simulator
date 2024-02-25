#include "network.h"

const int TRIANGULAR_NODE_NET_CNXS = 6;
const int TRIANGULAR_NODE_DUAL_CNXS = 6;

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
Network::Network() : maxNetCnxs(0), maxDualCnxs(0), dimensions(2), reciprocalDimensions(2), geometryCode("") {
    nodes.reserve(0);
}
/**
 * @brief Construct a network with a given number of nodes and max coordination
 * @param numNodes Number of nodes
 * @param maxCnxs Maximum number of connections (net and dual)
 */
Network::Network(const int &numNodes, const int &maxCnxs) {
    nodes.reserve(numNodes);
    for (int i = 0; i < numNodes; ++i) {
        Node node(i, maxCnxs, maxCnxs);
        nodes.push_back(node);
    }
}
/**
 * @brief Construct a triangular lattice
 * @param nNodes Number of nodes
 */
Network::Network(const int &nNodes) {
    initialiseTriangularLattice(roundedSqrt(nNodes));
    initialiseDescriptors(TRIANGULAR_NODE_NET_CNXS);
}

/**
 * @brief Initialise node and edge distribution
 * @param maxCnxs Maximum number of connections
 * @param logger Logger to log to
 */
void Network::initialiseDescriptors(const int &maxCnxs) {
    // Set sizes of vectors and matrices
    nodeDistribution = std::vector<int>(maxCnxs + 1);
    edgeDistribution = std::vector<std::vector<int>>(maxCnxs + 1);
    for (int i = 0; i < maxCnxs + 1; ++i)
        edgeDistribution[i] = std::vector<int>(maxCnxs + 1);

    // Count number of each node type and add to vector
    std::for_each(nodes.begin(), nodes.end(), [this](const Node &node) {
        ++nodeDistribution[node.netCnxs.size()];
    });

    // Double count number of each edge type and add to vector
    std::for_each(nodes.begin(), nodes.end(), [this](const Node &node) {
        int netCnxs_n_i = node.netCnxs.size();
        for (int j = 0; j < netCnxs_n_i; ++j) {
            int netCnxs_j = node.netCnxs[j];
            int netCnxs_n_j = nodes[netCnxs_j].netCnxs.size();
            ++edgeDistribution[netCnxs_n_i][netCnxs_n_j];
        }
    });
}

/**
 * @brief Construct a network from files
 * @param prefix Prefix of files to load
 * @param maxBaseCoordinationArg Maximum base coordination of nodes
 * @param maxDualCoordinationArg Maximum dual coordination of nodes
 * @param logger Logger to log to
 * @throw std::runtime_error if cannot open aux file, crds file, net file or dual file
 */
Network::Network(const std::string &pathPrefix, const int &maxBaseCoordinationArg,
                 const int &maxDualCoordinationArg, LoggerPtr logger) {
    logger->debug("Reading aux file {} ...", pathPrefix + "_aux.dat");

    // Initialise variables with aux file information
    std::istringstream ss("");
    std::ifstream auxFile(pathPrefix + "_aux.dat", std::ios::in);
    if (!auxFile.is_open()) {
        logger->critical("Aux file not found!");
        throw std::runtime_error("Aux file not found!");
    }

    int nNodes;
    std::string line;
    getline(auxFile, line);
    std::istringstream(line) >> nNodes;
    logger->debug("Number of nodes: {}", nNodes);
    getline(auxFile, line);
    ss.str(line);
    ss >> maxNetCnxs;
    ss >> maxDualCnxs;

    if (maxNetCnxs < maxBaseCoordinationArg)
        maxNetCnxs = maxBaseCoordinationArg;
    if (maxDualCnxs < maxDualCoordinationArg)
        maxDualCnxs = maxDualCoordinationArg;

    if (maxDualCnxs > 12) {
        maxNetCnxs += 20;
        maxDualCnxs += 20;
    }
    logger->debug("Max Net/Dual Connections: {} {}", maxNetCnxs, maxDualCnxs);
    getline(auxFile, line);
    std::istringstream(line) >> geometryCode;
    dimensions = std::vector<double>(2);
    reciprocalDimensions = std::vector<double>(2);
    getline(auxFile, line);
    ss.str(line);
    ss >> dimensions[0];
    ss >> dimensions[1];
    getline(auxFile, line);
    ss.str(line);
    ss >> reciprocalDimensions[0];
    ss >> reciprocalDimensions[1];
    auxFile.close();

    nodes.reserve(nNodes);
    for (int i = 0; i < nNodes; ++i) {
        Node node(i, maxNetCnxs, maxDualCnxs);
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
    logger->debug("Max net connections: {}", maxNetCnxs);
    logger->debug("Number of nodes: {}", nodes.size());
    initialiseDescriptors(maxNetCnxs);
}

// Initialise triangular lattice of periodic 6-coordinate nodes
void Network::initialiseTriangularLattice(const int &dim) {
    geometryCode = "2DE"; // 2D euclidean
    int dimSq = dim * dim;

    // make 6 coordinate nodes
    nodes.reserve(dimSq);
    for (int i = 0; i < nodes.size(); ++i) {
        Node node(i, TRIANGULAR_NODE_NET_CNXS, TRIANGULAR_NODE_DUAL_CNXS);
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
    Network dualNetwork(numberOfUniqueDualNodes, maxCnxs);

    addUnorderedDualConnections(dualNetwork);
    orderDualConnections(dualNetwork);
    addOrderedNetworkConnections(dualNetwork);
    centreRings(dualNetwork);

    // set remaining parameters
    dualNetwork.dimensions = dimensions;
    dualNetwork.reciprocalDimensions = reciprocalDimensions;
    dualNetwork.geometryCode = geometryCode;
    dualNetwork.initialiseDescriptors(maxCnxs);

    return dualNetwork;
}
/**
 * @brief Find the number of unique dual nodes in the entire network
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
        std::vector<int> orderedCnxs(maxDualCnxs);
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

void Network::centreRings(const Network &baseNetwork) {
    // Sync B coordinates
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

// Rescale coordinates and lattice dimensions
void Network::rescale(const double &scaleFactor) {
    vectorMultiply(dimensions, scaleFactor);
    vectorDivide(reciprocalDimensions, scaleFactor);
    std::for_each(nodes.begin(), nodes.end(), [&scaleFactor](Node &node) {
        vectorMultiply(node.crd, scaleFactor);
    });
}

// Get proportion of nodes of each size
std::vector<double> Network::getNodeDistribution() const {
    std::vector<double> nodeDistributionDouble(nodeDistribution.begin(), nodeDistribution.end());
    vectorNormalise(nodeDistributionDouble);
    return nodeDistributionDouble;
}

// Get proportion of edges with a node of each size at each end
std::vector<std::vector<double>> Network::getEdgeDistribution() const {
    std::vector<std::vector<double>> normalisedDist(edgeDistribution.size());
    double sum = 0.0;
    for (int i = 0; i < edgeDistribution.size(); ++i) {
        normalisedDist[i] = std::vector<double>(edgeDistribution[i].size());
        for (int j = 0; j < edgeDistribution[i].size(); ++j)
            normalisedDist[i][j] = edgeDistribution[i][j];
        sum += vectorSum(edgeDistribution[i]);
    }
    for (int i = 0; i < edgeDistribution.size(); ++i)
        vectorDivide(normalisedDist[i], sum);

    return normalisedDist;
}

// Calculate Aboav-Weaire fitting parameters
std::tuple<double, double, double> Network::getAboavWeaireParams() const {
    /* Aboav-Weaire's law: nm_n=<n>^2+mu+<n>(1-alpha)(n-<n>)
     * calculate mean node size, <n>^2
     * calculate mean node size about a node of size n
     * perform linear fit */

    // mean from node distribution
    double mean = 0.0;
    for (int i = 0; i < nodeDistribution.size(); ++i)
        mean += i * nodeDistribution[i];
    mean /= vectorSum(nodeDistribution);

    // find x,y only for sizes which are present
    std::vector<double> x(edgeDistribution.size());
    std::vector<double> y(edgeDistribution.size());
    for (int i = 0; i < edgeDistribution.size(); ++i) {
        const int num = vectorSum(edgeDistribution[i]);
        if (num > 0) {
            double mn = 0.0;
            for (int j = 0; j < edgeDistribution[i].size(); ++j)
                mn += j * edgeDistribution[i][j];
            mn /= num;
            x.push_back(mean * (i - mean));
            y.push_back(i * mn);
        }
    }

    // linear fit if more than one data point, return: alpha, mu, rsq
    if (x.size() > 1) {
        double gradient;
        double intercept;
        double rSquared;
        std::tie(gradient, intercept, rSquared) = vectorLinearRegression(x, y);
        return std::make_tuple(1.0 - gradient, intercept - mean * mean, rSquared);
    }
    return std::make_tuple(0.0, 0.0, 0.0);
}

// Calculate network assortativity through the degree correlation coefficient
double Network::getAssortativityOld() const {
    /* definitions:
     * 1) e_ij degree correlation matrix, prob of finding nodes with degree i,j at
     * end of random link 2) q_k prob of finding node with degree k at end of
     * random link, q_k=kp_k/<k> */
    std::vector<std::vector<double>> e = getEdgeDistribution();
    std::vector<double> q = getNodeDistribution();
    double mean = 0.0;
    for (int i = 0; i < q.size(); ++i)
        mean += i * q[i];
    for (int i = 0; i < q.size(); ++i)
        q[i] = i * q[i] / mean;

    /* degree correlation coefficient:
     * r=sum_jk jk(e_jk-q_j*q_k)/sig^2
     * sig^2 = sum_k k^2*q_k - (sum_k k*q_k)^2
     * bounded between -1 (perfectly disasssortative) and 1 (perfectly
     * assortative) with 0 as neutral */
    double r = 0.0;
    for (int j = 0; j < e.size(); ++j) {
        for (int k = 0; k < e.size(); ++k)
            r += j * k * (e[j][k] - q[j] * q[k]);
    }
    double sigSq;
    double a = 0.0;
    double b = 0.0; // dummy variables
    for (int k = 0; k < q.size(); ++k) {
        a += k * k * q[k];
        b += k * q[k];
    }
    sigSq = a - b * b;
    r /= sigSq;

    return r;
}

// Estimate alpha parameter from degree correlation coefficient
double Network::getAboavWeaireEstimate() {
    /* Can derive by substituting aw law into equation for r */
    std::vector<double> p = getNodeDistribution();
    std::vector<double> k(p.size());
    for (int i = 0; i < p.size(); ++i)
        k[i] = i;
    double n;
    double n2;
    double n3;
    double nSq;
    n = vectorSum(multiplyVectors(k, p));
    n2 = vectorSum(multiplyVectors(multiplyVectors(k, k), p));
    n3 = vectorSum(multiplyVectors(multiplyVectors(multiplyVectors(k, k), k), p));
    nSq = n * n;
    double alpha;
    double r = getAssortativity();
    double mu = n2 - nSq;
    alpha = (-r * (n * n3 - n2 * n2) - mu * mu) / (nSq * mu);

    return alpha;
}

// Calculate entropy of node and edge distribution
std::vector<double> Network::getEntropyOld() const {
    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    std::vector<double> p = getNodeDistribution();
    std::vector<double> q = p;
    double mean = 0.0;
    for (int i = 0; i < q.size(); ++i)
        mean += i * q[i];
    for (int i = 0; i < q.size(); ++i)
        q[i] = i * q[i] / mean;
    std::vector<std::vector<double>> e = getEdgeDistribution();

    std::for_each(p.begin(), p.end(), [&s0](const double &i) {
        if (i > 0.0)
            s0 -= i * log(i);
    });

    for (int i = 0; i < e.size(); ++i) {
        for (int j = 0; j < e[i].size(); ++j) {
            if (e[i][j] > 0.0) {
                s1 -= e[i][j] * log(e[i][j]);
                s2 += e[i][j] * log(e[i][j] / (q[i] * q[j]));
            }
        }
    }

    return std::vector<double>{s0, s1, s2};
}

double Network::getAboavWeaire() const {
    return 0.0;
}
double Network::getEntropy() const {
    return 0.0;
}

double Network::getAssortativity() const {
    return 0.0;
}

/**
 * @brief Get the number of nodes with a given number of connections
 * @param minRingSize Minimum number of connections
 * @param maxRingSize Maximum number of connections
 * @return Vector of the number of nodes with a given number of connections
 */
std::vector<int> Network::getCoordinations(const int &minCoordination, const int &maxCoordination) const {
    std::vector<int> ringSizes(maxCoordination - minCoordination + 1, 0);
    for (int i = 0; i < nodes.size(); ++i) {
        try {
            ringSizes[nodes[i].netCnxs.size() - minCoordination]++;
        } catch (std::out_of_range &e) {
            throw std::runtime_error("Node " + std::to_string(i) + " has " + std::to_string(nodes[i].netCnxs.size()) +
                                     " connections, which is outside the range " + std::to_string(minCoordination) +
                                     " to " + std::to_string(maxCoordination));
        }
    }
    return ringSizes;
}

// Write network in format which can be loaded
void Network::write(const std::string &prefix) {
    // auxilary information
    std::ofstream auxFile(prefix + "_aux.dat", std::ios::in | std::ios::trunc);
    auxFile << std::fixed << std::showpoint << std::setprecision(1);
    auxFile << std::setw(10) << std::left << nodes.size() << std::endl;
    auxFile << std::setw(10) << std::left << getMaxCnxs()
            << std::setw(10) << std::left << getMaxDualCnxs() << std::endl;
    auxFile << std::setw(10) << std::left << geometryCode << std::endl;
    auxFile << std::fixed << std::showpoint << std::setprecision(6);
    for (int i = 0; i < dimensions.size(); ++i)
        auxFile << std::setw(20) << std::left << dimensions[i];
    auxFile << std::endl;
    for (int i = 0; i < reciprocalDimensions.size(); ++i)
        auxFile << std::setw(20) << std::left << reciprocalDimensions[i];
    auxFile << std::endl;
    auxFile.close();

    // coordinates
    std::ofstream crdFile(prefix + "_crds.dat", std::ios::in | std::ios::trunc);
    crdFile << std::fixed << std::showpoint << std::setprecision(6);
    for (int i = 0; i < nodes.size(); ++i) {
        for (int j = 0; j < nodes[i].crd.size(); ++j) {
            crdFile << std::setw(20) << std::left << nodes[i].crd[j];
        }
        crdFile << std::endl;
    }
    crdFile.close();

    // network connections
    std::ofstream netFile(prefix + "_net.dat", std::ios::in | std::ios::trunc);
    netFile << std::fixed << std::showpoint << std::setprecision(1);
    for (int i = 0; i < nodes.size(); ++i) {
        for (int j = 0; j < nodes[i].netCnxs.size(); ++j) {
            netFile << std::setw(20) << std::left << nodes[i].netCnxs[j];
        }
        netFile << std::endl;
    }
    netFile.close();

    // dual connections
    std::ofstream dualFile(prefix + "_dual.dat", std::ios::in | std::ios::trunc);
    dualFile << std::fixed << std::showpoint << std::setprecision(1);
    for (int i = 0; i < nodes.size(); ++i) {
        for (int j = 0; j < nodes[i].dualCnxs.size(); ++j) {
            dualFile << std::setw(20) << std::left << nodes[i].dualCnxs[j];
        }
        dualFile << std::endl;
    }
    dualFile.close();
}

int Network::getMaxCnxs() {
    auto maxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.netCnxs.size() < b.netCnxs.size();
    });
    return maxNode->netCnxs.size();
}

int Network::getMinCnxs() {
    auto minNode = std::min_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.netCnxs.size() < b.netCnxs.size();
    });
    return minNode->netCnxs.size();
}

int Network::getMaxDualCnxs() {
    auto maxNode = std::max_element(nodes.begin(), nodes.end(), [](const Node &a, const Node &b) {
        return a.dualCnxs.size() < b.dualCnxs.size();
    });
    return maxNode->dualCnxs.size();
}

int Network::getMinDualCnxs() {
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
