#include "network.h"
#include <cmath>

const int TRIANGULAR_NODE_NET_CNXS = 6;
const int TRIANGULAR_NODE_DUAL_CNXS = 6;

/**
 * @brief Calculate the rounded square root of a number
 * @param num Number to calculate the rounded square root of
 * @return Rounded square root of the number
 */
inline int roundedSqrt(int num) {
    return std::round(std::sqrt(num));
}

// Default constructor
Network::Network() : nodes(VecR<Node>(0, 1)) {}

// Construct with number of nodes and connectivity limits
Network::Network(int nNodes, int maxCnxs) : nodes(VecR<Node>(0, nNodes)) {
    for (int i = 0; i < nNodes; ++i) {
        Node node(i, maxCnxs, maxCnxs, 0);
        nodes.addValue(node);
    }
}
/**
 * @brief Construct a triangular lattice
 * @param nNodes Number of nodes
 */
Network::Network(int nNodes) {
    initialiseTriangularLattice(roundedSqrt(nNodes));
    initialiseDescriptors(TRIANGULAR_NODE_NET_CNXS);
}

/**
 * @brief Initialise node and edge distribution
 * @param maxCnxs Maximum number of connections
 * @param logger Logger to log to
 */
void Network::initialiseDescriptors(int maxCnxs) {
    // Set sizes of vectors and matrices
    nodeDistribution = VecF<int>(maxCnxs + 1);
    edgeDistribution = VecF<VecF<int>>(maxCnxs + 1);
    for (int i = 0; i < maxCnxs + 1; ++i)
        edgeDistribution[i] = VecF<int>(maxCnxs + 1);

    // Count number of each node type and add to vector
    for (int i = 0; i < nodes.n; ++i)
        ++nodeDistribution[nodes[i].netCnxs.n];

    // Double count number of each edge type and add to vector
    for (int i = 0; i < nodes.n; ++i) {
        int netCnxs_n_i = nodes[i].netCnxs.n;
        for (int j = 0; j < netCnxs_n_i; ++j) {
            int netCnxs_j = nodes[i].netCnxs[j];
            int netCnxs_n_j = nodes[netCnxs_j].netCnxs.n;
            ++edgeDistribution[netCnxs_n_i][netCnxs_n_j];
        }
    }
}

/**
 * @brief Construct a network from files
 * @param prefix Prefix of files to load
 * @param maxBaseCoordinationArg Maximum base coordination of nodes
 * @param maxDualCoordinationArg Maximum dual coordination of nodes
 * @param logger Logger to log to
 */
Network::Network(const std::string &prefix, int maxBaseCoordinationArg,
                 int maxDualCoordinationArg, LoggerPtr logger) {
    logger->info("Reading aux file {} ...", prefix + "_aux.dat");

    // Initialise variables with aux file information
    std::string line;
    std::istringstream ss("");
    std::ifstream auxFile(prefix + "_aux.dat", std::ios::in);
    if (!auxFile.is_open()) {
        logger->critical("Aux file not found!");
        throw std::runtime_error("Aux file not found!");
    }

    int nNodes;
    getline(auxFile, line);
    std::istringstream(line) >> nNodes;
    logger->info("Number of nodes: {}", nNodes);
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
    logger->info("Max Net/Dual Connections: {} {}", maxNetCnxs, maxDualCnxs);
    getline(auxFile, line);
    std::istringstream(line) >> geometryCode;
    nodes = VecR<Node>(0, nNodes);
    dimensions = VecF<double>(2);
    rpb = VecF<double>(2);
    getline(auxFile, line);
    ss.str(line);
    ss >> dimensions[0];
    ss >> dimensions[1];
    getline(auxFile, line);
    ss.str(line);
    ss >> rpb[0];
    ss >> rpb[1];
    for (int i = 0; i < nodes.nMax; ++i) {
        Node node(i, maxNetCnxs, maxDualCnxs, 0);
        nodes.addValue(node);
    }
    auxFile.close();

    logger->info("Reading crds file {} ...", prefix + "_crds.dat");
    std::ifstream crdFile(prefix + "_crds.dat", std::ios::in);
    if (!crdFile.is_open()) {
        logger->critical("crds.dat file not found!");
        throw std::runtime_error("crds.dat file not found!");
    }
    VecF<double> crd(2);
    for (int i = 0; i < nodes.n; ++i) {
        getline(crdFile, line);
        ss.str(line);
        ss >> crd[0];
        ss >> crd[1];
        nodes[i].crd = crd;
    }
    crdFile.close();

    // Read network connections
    logger->info("Reading net file {} ...", prefix + "_net.dat");
    std::ifstream netFile(prefix + "_net.dat", std::ios::in);
    if (!netFile.is_open()) {
        logger->critical("net.dat file not found!");
        throw std::runtime_error("net.dat file not found!");
    }
    int cnx;
    for (int i = 0; i < nodes.n; ++i) {
        getline(netFile, line);
        std::istringstream ss(line);
        while (ss >> cnx) {
            nodes[i].netCnxs.addValue(cnx);
        }
    }
    netFile.close();

    // Read dual connections
    logger->info("Reading dual file {} ...", prefix + "_dual.dat");
    std::ifstream dualFile(prefix + "_dual.dat", std::ios::in);
    if (!dualFile.is_open()) {
        logger->critical("dual.dat file not found!");
        throw std::runtime_error("dual.dat file not found!");
    }
    for (int i = 0; i < nodes.n; ++i) {
        getline(dualFile, line);
        std::istringstream ss(line);
        while (ss >> cnx)
            nodes[i].dualCnxs.addValue(cnx);
    }
    dualFile.close();

    // Set up descriptors
    logger->info("Max net connections: {}", maxNetCnxs);
    logger->info("Number of nodes: {}", nodes.n);
    initialiseDescriptors(maxNetCnxs);
}

// Initialise triangular lattice of periodic 6-coordinate nodes
void Network::initialiseTriangularLattice(int dim) {
    geometryCode = "2DE"; // 2D euclidean
    int dimSq = dim * dim;
    nodes = VecR<Node>(0, dimSq);

    // make 6 coordinate nodes
    for (int i = 0; i < nodes.nMax; ++i) {
        Node node(i, TRIANGULAR_NODE_NET_CNXS, TRIANGULAR_NODE_DUAL_CNXS, 0);
        nodes.addValue(node);
    }

    // assign coordinates in layers, with unit bond lengths
    dimensions = VecF<double>(2);
    rpb = VecF<double>(2);
    dimensions[0] = dim;
    dimensions[1] = dim * sqrt(3) * 0.5;
    rpb[0] = 1.0 / dimensions[0];
    rpb[1] = 1.0 / dimensions[1];
    assignCoordinates(dim);
    makeConnections(dim, dimSq);
    makeDualConnections(dim, dimSq);
}

void Network::assignCoordinates(int dim) {
    VecF<double> c(2);
    double dy = sqrt(3.0) * 0.5;
    for (int y = 0; y < dim; ++y) {
        c[1] = 0.5 * dy + y * dy;
        for (int x = 0; x < dim; ++x) {
            c[0] = 0.5 * (y % 2) + x;
            nodes[x + y * dim].crd = c;
        }
    }
}

void Network::makeConnections(int dim, int dimSq) {
    int id = 0;
    int cnx;
    for (int y = 0; y < dim; ++y) {
        for (int x = 0; x < dim; ++x) {
            cnx = y * dim + (id + dim - 1) % dim;
            nodes[id].netCnxs.addValue(cnx);
            addConnectionsBasedOnParity(y, id, cnx, dim, dimSq);
            addMoreConnectionsBasedOnParity(y, id, cnx, dim, dimSq);
            ++id;
        }
    }
}

void Network::addConnectionsBasedOnParity(int y, int id, int cnx, int dim,
                                          int dimSq) {
    if (y % 2 == 0) {
        cnx = (cnx + dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
        cnx = (id + dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
    } else {
        cnx = (id + dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
        cnx = ((y + 1) * dim + (id + 1) % dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
    }
    cnx = y * dim + (id + 1) % dim;
    nodes[id].netCnxs.addValue(cnx);
}

void Network::addMoreConnectionsBasedOnParity(int y, int id, int cnx, int dim,
                                              int dimSq) {
    if (y % 2 == 0) {
        cnx = (id + dimSq - dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
        cnx = (dimSq + (y - 1) * dim + (id + dim - 1) % dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
    } else {
        cnx = (dimSq + (y - 1) * dim + (id + dim + 1) % dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
        cnx = (dimSq + (y - 1) * dim + (id + dim) % dim) % dimSq;
        nodes[id].netCnxs.addValue(cnx);
    }
}

void Network::makeDualConnections(int dim, int dimSq) {
    int id = 0;
    int dimSq2 = 2 * dimSq;
    for (int y = 0; y < dim; ++y) {
        for (int x = 0; x < dim; ++x) {
            addDualConnections(y, id, dim, dimSq2);
            ++id;
        }
    }
}

void Network::addDualConnections(int y, int id, int dim, int dimSq2) {
    int cnx = (2 * y * dim + (id + dim - 1) % dim);
    nodes[id].dualCnxs.addValue(cnx);
    cnx = ((2 * y + 1) * dim + (id) % dim);
    nodes[id].dualCnxs.addValue(cnx);
    cnx = (2 * y * dim + id % dim);
    nodes[id].dualCnxs.addValue(cnx);

    if (y % 2 == 0) {
        cnx = (2 * y * dim + id % dim - dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.addValue(cnx);
        cnx = (2 * (y - 1) * dim + (id + dim - 1) % dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.addValue(cnx);
        cnx = (2 * (y - 1) * dim + (id + dim - 1) % dim + dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.addValue(cnx);
    } else {
        cnx = (2 * y * dim + (id + 1) % dim - dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.addValue(cnx);
        cnx = (2 * (y - 1) * dim + id % dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.addValue(cnx);
        cnx = (2 * (y - 1) * dim + id % dim + dim + dimSq2) % dimSq2;
        nodes[id].dualCnxs.addValue(cnx);
    }
}

Network Network::constructDual(int maxCnxs) {
    int numberOfUniqueDualNodes = findNumberOfUniqueDualNodes();
    Network dualNetwork(numberOfUniqueDualNodes, maxCnxs);

    addUnorderedDualConnections(dualNetwork);
    orderDualConnections(dualNetwork);
    addOrderedNetworkConnections(dualNetwork);
    setCoordinatesAtCentreOfDualConnections(dualNetwork);

    // set remaining parameters
    dualNetwork.dimensions = dimensions;
    dualNetwork.rpb = rpb;
    dualNetwork.geometryCode = geometryCode;
    dualNetwork.initialiseDescriptors(maxCnxs);

    return dualNetwork;
}

int Network::findNumberOfUniqueDualNodes() {
    int numberOfUniqueDualNodes = -1;
    for (int i = 0; i < nodes.n; ++i) {
        for (int j = 0; j < nodes[i].dualCnxs.n; ++j) {
            if (nodes[i].dualCnxs[j] > numberOfUniqueDualNodes)
                numberOfUniqueDualNodes = nodes[i].dualCnxs[j];
        }
    }
    return numberOfUniqueDualNodes + 1;
}

void Network::addUnorderedDualConnections(Network &dualNetwork) {
    for (int i = 0; i < nodes.n; ++i) {
        for (int j = 0; j < nodes[i].dualCnxs.n; ++j) {
            int dualConnectionId = nodes[i].dualCnxs[j];
            dualNetwork.nodes[dualConnectionId].dualCnxs.addValue(i);
        }
    }
}

void Network::orderDualConnections(Network &dualNetwork) {
    for (int i = 0; i < dualNetwork.nodes.n; ++i) {
        VecR<int> unordered = dualNetwork.nodes[i].dualCnxs;
        VecR<int> ordered(0, unordered.nMax);
        ordered.addValue(unordered[0]);
        for (int j = 1; j < unordered.n; ++j) {
            VecR<int> common =
                vCommonValues(nodes[ordered[j - 1]].netCnxs, unordered);
            if (!vContains(ordered, common[0]))
                ordered.addValue(common[0]);
            else
                ordered.addValue(common[1]);
        }
        dualNetwork.nodes[i].dualCnxs = ordered;
    }
}

void Network::addOrderedNetworkConnections(Network &dualNetwork) {
    for (int i = 0; i < dualNetwork.nodes.n; ++i) {
        VecR<int> dualCnxs = dualNetwork.nodes[i].dualCnxs;
        VecR<int> common;
        for (int j = 0; j < dualCnxs.n; ++j) {
            int k = (j + 1) % dualCnxs.n;
            common = vCommonValues(nodes[dualCnxs[j]].dualCnxs,
                                   nodes[dualCnxs[k]].dualCnxs);
            common.delValue(i);
            dualNetwork.nodes[i].netCnxs.addValue(common[0]);
        }
    }
}

void Network::setCoordinatesAtCentreOfDualConnections(Network &dualNetwork) {
    for (int i = 0; i < dualNetwork.nodes.n; ++i) {
        VecF<double> x(dualNetwork.nodes[i].dualCnxs.n);
        VecF<double> y(dualNetwork.nodes[i].dualCnxs.n);
        for (int j = 0; j < dualNetwork.nodes[i].dualCnxs.n; ++j) {
            x[j] = nodes[dualNetwork.nodes[i].dualCnxs[j]].crd[0];
            y[j] = nodes[dualNetwork.nodes[i].dualCnxs[j]].crd[1];
        }
        VecF<double> origin(2);
        origin[0] = x[0];
        origin[1] = y[0];
        x -= origin[0];
        y -= origin[1];
        for (int j = 0; j < x.n; ++j)
            x[j] -= dimensions[0] * nearbyint(x[j] * rpb[0]);
        for (int j = 0; j < y.n; ++j)
            y[j] -= dimensions[1] * nearbyint(y[j] * rpb[1]);
        VecF<double> c(2);
        c[0] = origin[0] + vMean(x);
        c[1] = origin[1] + vMean(y);
        dualNetwork.nodes[i].crd = c;
    }
}

// Generate auxilary connections
void Network::gen2ndOrderConnections(Network dualNetwork) {

    // generate second order network connections (share single point in dual)
    // increase size of aux containers
    for (int i = 0; i < nodes.n; ++i)
        nodes[i].auxCnxs = VecR<int>(0, nodes[i].dualCnxs.nMax);

    // generate unordered second order connections
    for (int i = 0; i < dualNetwork.nodes.n; ++i) {
        VecR<int> dualCnxs = dualNetwork.nodes[i].dualCnxs;
        for (int j = 0; j < dualCnxs.n - 1; ++j) {
            int id0 = dualCnxs[j];
            for (int k = j + 1; k < dualCnxs.n; ++k) {
                int id1 = dualCnxs[k];
                if (!vContains(nodes[id0].netCnxs, id1) &&
                    !vContains(nodes[id0].auxCnxs, id1)) {
                    nodes[id0].auxCnxs.addValue(id1);
                    nodes[id1].auxCnxs.addValue(id0);
                }
            }
        }
    }
}

// Rescale coordinates and lattice dimensions
void Network::rescale(double scaleFactor) {
    dimensions *= scaleFactor;
    rpb /= scaleFactor;
    for (int i = 0; i < nodes.n; ++i)
        nodes[i].crd *= scaleFactor;
}

// Find local region of lattice, nodes in a given range of central nodes
void Network::findLocalRegion(int a, int b, int extent, VecR<int> &local,
                              VecR<int> &fixedInner, VecR<int> &fixedOuter) {

    VecR<int> localNodes(0, 10000);
    VecR<int> fixedInnerNodes(0, 10000);
    VecR<int> fixedOuterNodes(0, 10000);
    VecR<int> layer0(0, 2);
    VecR<int> layer1(0, 100);
    VecR<int> layer2;
    VecR<int> common;
    layer0.addValue(a);
    layer0.addValue(b);
    for (int i = 0; i < layer0.n; ++i)
        localNodes.addValue(layer0[i]);
    for (int i = 0; i < layer0.n; ++i) {
        for (int j = 0; j < nodes[layer0[i]].netCnxs.n; ++j) {
            layer1.addValue(nodes[layer0[i]].netCnxs[j]);
        }
    }
    layer1 = vUnique(layer1);
    for (int i = 0; i < layer0.n; ++i)
        layer1.delValue(layer0[i]);
    for (int i = 0; i < layer1.n; ++i)
        localNodes.addValue(layer1[i]);

    for (int i = 0; i < extent + 1; ++i) {
        layer2 = VecR<int>(0, 1000);
        for (int j = 0; j < layer1.n; ++j) {
            for (int k = 0; k < nodes[layer1[j]].netCnxs.n; ++k) {
                layer2.addValue(nodes[layer1[j]].netCnxs[k]);
            }
        }
        if (layer2.n == 0)
            break;
        layer2 = vUnique(layer2);
        common = vCommonValues(layer2, layer1);
        for (int j = 0; j < common.n; ++j)
            layer2.delValue(common[j]);
        common = vCommonValues(layer2, layer0);
        for (int j = 0; j < common.n; ++j)
            layer2.delValue(common[j]);
        if (i < extent - 1) {
            for (int j = 0; j < layer2.n; ++j)
                localNodes.addValue(layer2[j]);
        } else if (i == extent - 1) {
            for (int j = 0; j < layer2.n; ++j)
                fixedInnerNodes.addValue(layer2[j]);
        } else if (i == extent) {
            for (int j = 0; j < layer2.n; ++j)
                fixedOuterNodes.addValue(layer2[j]);
        }
        layer0 = layer1;
        layer1 = layer2;
    }

    local = VecR<int>(localNodes.n);
    fixedInner = VecR<int>(fixedInnerNodes.n);
    fixedOuter = VecR<int>(fixedOuterNodes.n);
    for (int i = 0; i < local.n; ++i)
        local[i] = localNodes[i];
    for (int i = 0; i < fixedInner.n; ++i)
        fixedInner[i] = fixedInnerNodes[i];
    for (int i = 0; i < fixedOuter.n; ++i)
        fixedOuter[i] = fixedOuterNodes[i];
}

// Get proportion of nodes of each size
VecF<double> Network::getNodeDistribution() {

    VecF<double> normalisedDist(nodeDistribution.n);
    for (int i = 0; i < nodeDistribution.n; ++i)
        normalisedDist[i] = nodeDistribution[i];
    normalisedDist /= vSum(normalisedDist);
    return normalisedDist;
}

// Get proportion of edges with a node of each size at each end
VecF<VecF<double>> Network::getEdgeDistribution() {

    VecF<VecF<double>> normalisedDist(edgeDistribution.n);
    double sum = 0.0;
    for (int i = 0; i < edgeDistribution.n; ++i) {
        normalisedDist[i] = VecF<double>(edgeDistribution[i].n);
        for (int j = 0; j < edgeDistribution[i].n; ++j)
            normalisedDist[i][j] = edgeDistribution[i][j];
        sum += vSum(edgeDistribution[i]);
    }
    for (int i = 0; i < edgeDistribution.n; ++i)
        normalisedDist[i] /= sum;

    return normalisedDist;
}

// Calculate Aboav-Weaire fitting parameters
VecF<double> Network::aboavWeaireParams() {

    /* Aboav-Weaire's law: nm_n=<n>^2+mu+<n>(1-alpha)(n-<n>)
     * calculate mean node size, <n>^2
     * calculate mean node size about a node of size n
     * perform linear fit */

    // mean from node distribution
    double mean = 0.0;
    for (int i = 0; i < nodeDistribution.n; ++i)
        mean += i * nodeDistribution[i];
    mean /= vSum(nodeDistribution);

    // find x,y only for sizes which are present
    VecR<double> x(0, edgeDistribution.n);
    VecR<double> y(0, edgeDistribution.n);
    for (int i = 0; i < edgeDistribution.n; ++i) {
        int num = vSum(edgeDistribution[i]);
        if (num > 0) {
            double mn = 0.0;
            for (int j = 0; j < edgeDistribution[i].n; ++j)
                mn += j * edgeDistribution[i][j];
            mn /= num;
            x.addValue(mean * (i - mean));
            y.addValue(i * mn);
        }
    }

    // linear fit if more than one data point, return: alpha, mu, rsq
    VecF<double> aw(3);
    if (x.n > 1) {
        VecR<double> fit = vLinearRegression(x, y);
        aw[0] = 1.0 - fit[0];
        aw[1] = fit[1] - mean * mean;
        aw[2] = fit[2];
    }

    return aw;
}

// Calculate network assortativity through the degree correlation coefficient
double Network::assortativity() {

    /* definitions:
     * 1) e_ij degree correlation matrix, prob of finding nodes with degree i,j at
     * end of random link 2) q_k prob of finding node with degree k at end of
     * random link, q_k=kp_k/<k> */
    VecF<VecF<double>> e = getEdgeDistribution();
    VecF<double> q = getNodeDistribution();
    double mean = 0.0;
    for (int i = 0; i < q.n; ++i)
        mean += i * q[i];
    for (int i = 0; i < q.n; ++i)
        q[i] = i * q[i] / mean;

    /* degree correlation coefficient:
     * r=sum_jk jk(e_jk-q_j*q_k)/sig^2
     * sig^2 = sum_k k^2*q_k - (sum_k k*q_k)^2
     * bounded between -1 (perfectly disasssortative) and 1 (perfectly
     * assortative) with 0 as neutral */
    double r = 0.0;
    for (int j = 0; j < e.n; ++j) {
        for (int k = 0; k < e.n; ++k)
            r += j * k * (e[j][k] - q[j] * q[k]);
    }
    double sigSq;
    double a = 0.0;
    double b = 0.0; // dummy variables
    for (int k = 0; k < q.n; ++k) {
        a += k * k * q[k];
        b += k * q[k];
    }
    sigSq = a - b * b;
    r /= sigSq;

    return r;
}

// Estimate alpha parameter from degree correlation coefficient
double Network::aboavWeaireEstimate() {

    /* Can derive by substituting aw law into equation for r */
    VecF<double> p = getNodeDistribution();
    VecF<double> k(p.n);
    for (int i = 0; i < p.n; ++i)
        k[i] = i;
    double n;
    double n2;
    double n3;
    double nSq;
    n = vSum(k * p);
    n2 = vSum(k * k * p);
    n3 = vSum(k * k * k * p);
    nSq = n * n;
    double alpha;
    double r = assortativity();
    double mu = n2 - nSq;
    alpha = (-r * (n * n3 - n2 * n2) - mu * mu) / (nSq * mu);

    return alpha;
}

// Calculate entropy of node and edge distribution
VecF<double> Network::entropy() {

    double s0 = 0.0;
    double s1 = 0.0;
    double s2 = 0.0;
    VecF<double> p = getNodeDistribution();
    VecF<double> q = p;
    double mean = 0.0;
    for (int i = 0; i < q.n; ++i)
        mean += i * q[i];
    for (int i = 0; i < q.n; ++i)
        q[i] = i * q[i] / mean;
    VecF<VecF<double>> e = getEdgeDistribution();

    for (int i = 0; i < p.n; ++i) {
        if (p[i] > 0.0)
            s0 -= p[i] * log(p[i]);
    }

    for (int i = 0; i < e.n; ++i) {
        for (int j = 0; j < e[i].n; ++j) {
            if (e[i][j] > 0.0) {
                s1 -= e[i][j] * log(e[i][j]);
                s2 += e[i][j] * log(e[i][j] / (q[i] * q[j]));
            }
        }
    }

    VecF<double> entropy(3);
    entropy[0] = s0;
    entropy[1] = s1;
    entropy[2] = s2;

    return entropy;
}

// Get cluster statistics for given node coordination
VecF<int> Network::maxClusters(int minCnd, int maxCnd, int minInnerCnxs,
                               int minOuterCnxs) {

    // Loop over each coordination size
    VecF<int> maxClusters((maxCnd - minCnd) + 1);
    for (int cnd = minCnd, cndIndex = 0; cnd <= maxCnd; ++cnd, ++cndIndex) {

        // Identify nodes with the required coordination and similarly coordinated
        // neighbours
        VecF<int> innerNodes(nodes.n);
        VecF<int> outerNodes(nodes.n);
        for (int i = 0; i < nodes.n; ++i) {
            if (nodes[i].netCnxs.n == cnd) {
                int nCnxs = 0;
                for (int j = 0; j < cnd; ++j)
                    if (nodes[nodes[i].netCnxs[j]].netCnxs.n == cnd)
                        nCnxs += 1;
                if (nCnxs >= minInnerCnxs) {
                    innerNodes[i] = 1;
                    outerNodes[i] = 0;
                } else if (nCnxs >= minOuterCnxs) {
                    innerNodes[i] = 0;
                    outerNodes[i] = 1;
                } else {
                    innerNodes[i] = 0;
                    outerNodes[i] = 0;
                }
            } else {
                innerNodes[i] = 0;
                outerNodes[i] = 0;
            }
        }

        // Find largest cluster
        int maxClstSize = 0;
        for (int i = 0; i < nodes.n; ++i) {
            if (innerNodes[i] == 1) {
                VecR<int> clst(0, nodes.n);
                VecR<int> search0(0, nodes.n);
                VecR<int> search1(0, nodes.n);
                VecR<int> search2(0, nodes.n);
                clst.addValue(i);
                search0.addValue(i);
                for (;;) {
                    for (int j = 0; j < search0.n; ++j) {
                        for (int k = 0; k < cnd; ++k) {
                            int id = nodes[search0[j]].netCnxs[k];
                            if (innerNodes[id])
                                search1.addValue(id);
                            else if (outerNodes[id])
                                search2.addValue(id);
                        }
                    }
                    search1 = vUnique(search1);
                    VecR<int> delValues = vCommonValues(clst, search1);
                    for (int j = 0; j < delValues.n; ++j)
                        search1.delValue(delValues[j]);
                    delValues = vCommonValues(clst, search2);
                    for (int j = 0; j < delValues.n; ++j)
                        search2.delValue(delValues[j]);
                    for (int j = 0; j < search1.n; ++j)
                        clst.addValue(search1[j]);
                    for (int j = 0; j < search2.n; ++j)
                        clst.addValue(search2[j]);
                    search0 = search1;
                    search1 = VecR<int>(0, nodes.n);
                    search2 = VecR<int>(0, nodes.n);
                    if (search0.n == 0)
                        break;
                }
                for (int j = 0; j < clst.n; ++j) {
                    innerNodes[clst[j]] = 0;
                    outerNodes[clst[j]] = 0;
                }
                if (clst.n > maxClstSize)
                    maxClstSize = clst.n;
            }
        }

        // Add to results vector
        maxClusters[cndIndex] = maxClstSize;
    }

    return maxClusters;
}

// Get cluster statistics for given node coordination
double Network::maxCluster(int nodeCnd) {

    // Identify nodes with the required coordination
    VecF<int> activeNodes(nodes.n);
    for (int i = 0; i < nodes.n; ++i) {
        if (nodes[i].netCnxs.n == nodeCnd)
            activeNodes[i] = 1;
        else
            activeNodes[i] = 0;
    }

    // Find largest cluster
    int maxClstSize = 0;
    for (int i = 0; i < nodes.n; ++i) {
        if (activeNodes[i] == 1) {
            VecR<int> clst(0, nodes.n);
            VecR<int> search0(0, nodes.n);
            VecR<int> search1(0, nodes.n);
            clst.addValue(i);
            search0.addValue(i);
            for (;;) {
                for (int j = 0; j < search0.n; ++j) {
                    for (int k = 0; k < nodeCnd; ++k) {
                        int id = nodes[search0[j]].netCnxs[k];
                        if (nodes[id].netCnxs.n == nodeCnd)
                            search1.addValue(id);
                    }
                }
                if (search1.n == 0)
                    break;
                search1 = vUnique(search1);
                VecR<int> delValues = vCommonValues(clst, search1);
                for (int j = 0; j < delValues.n; ++j)
                    search1.delValue(delValues[j]);
                search0 = search1;
                search1 = VecR<int>(0, nodes.n);
                for (int j = 0; j < search0.n; ++j)
                    clst.addValue(search0[j]);
                if (search0.n == 0)
                    break;
            }
            for (int j = 0; j < clst.n; ++j)
                activeNodes[clst[j]] = 0;
            if (clst.n > maxClstSize)
                maxClstSize = clst.n;
        }
    }

    return maxClstSize;
}

// Write network in format which can be loaded
void Network::write(const std::string &prefix) {
    // auxilary information
    std::ofstream auxFile(prefix + "_aux.dat", std::ios::in | std::ios::trunc);
    auxFile << std::fixed << std::showpoint << std::setprecision(1);
    auxFile << std::setw(10) << std::left << nodes.n << std::endl;
    auxFile << std::setw(10) << std::left << nodes[0].netCnxs.nMax
            << std::setw(10) << std::left << nodes[0].dualCnxs.nMax << std::endl;
    auxFile << std::setw(10) << std::left << geometryCode << std::endl;
    auxFile << std::fixed << std::showpoint << std::setprecision(6);
    for (int i = 0; i < dimensions.n; ++i)
        auxFile << std::setw(20) << std::left << dimensions[i];
    auxFile << std::endl;
    for (int i = 0; i < rpb.n; ++i)
        auxFile << std::setw(20) << std::left << rpb[i];
    auxFile << std::endl;
    auxFile.close();

    // coordinates
    std::ofstream crdFile(prefix + "_crds.dat", std::ios::in | std::ios::trunc);
    crdFile << std::fixed << std::showpoint << std::setprecision(6);
    for (int i = 0; i < nodes.n; ++i) {
        for (int j = 0; j < nodes[i].crd.n; ++j) {
            crdFile << std::setw(20) << std::left << nodes[i].crd[j];
        }
        crdFile << std::endl;
    }
    crdFile.close();

    // network connections
    std::ofstream netFile(prefix + "_net.dat", std::ios::in | std::ios::trunc);
    netFile << std::fixed << std::showpoint << std::setprecision(1);
    for (int i = 0; i < nodes.n; ++i) {
        for (int j = 0; j < nodes[i].netCnxs.n; ++j) {
            netFile << std::setw(20) << std::left << nodes[i].netCnxs[j];
        }
        netFile << std::endl;
    }
    netFile.close();

    // dual connections
    std::ofstream dualFile(prefix + "_dual.dat", std::ios::in | std::ios::trunc);
    dualFile << std::fixed << std::showpoint << std::setprecision(1);
    for (int i = 0; i < nodes.n; ++i) {
        for (int j = 0; j < nodes[i].dualCnxs.n; ++j) {
            dualFile << std::setw(20) << std::left << nodes[i].dualCnxs[j];
        }
        dualFile << std::endl;
    }
    dualFile.close();
}

// Write xyz file format of network
void Network::writeXYZ(const std::string &prefix, const std::string &element) {
    std::ofstream xyzFile(prefix + ".xyz", std::ios::in | std::ios::trunc);
    if (geometryCode == "2DE") {
        xyzFile << nodes.n << std::endl;
        xyzFile << "# Element X Y Z  " << std::endl;
        xyzFile << std::fixed << std::showpoint << std::setprecision(6);
        for (int i = 0; i < nodes.n; ++i) {
            xyzFile << std::setw(10) << std::left << element + std::to_string(i);
            xyzFile << std::setw(20) << std::left << nodes[i].crd[0];
            xyzFile << std::setw(20) << std::left << nodes[i].crd[1];
            xyzFile << std::setw(20) << std::left << "0.0" << std::endl;
        }
    } else if (geometryCode == "2DS") {
        xyzFile << nodes.n << std::endl;
        xyzFile << "# Element X Y Z  " << std::endl;
        xyzFile << std::fixed << std::showpoint << std::setprecision(6);
        for (int i = 0; i < nodes.n; ++i) {
            xyzFile << std::setw(10) << std::left << element + std::to_string(i);
            xyzFile << std::setw(20) << std::left << nodes[i].crd[0];
            xyzFile << std::setw(20) << std::left << nodes[i].crd[1];
            xyzFile << std::setw(20) << std::left << nodes[i].crd[2] << std::endl;
        }
    } else if (geometryCode == "2DEtr") {
        int nSi;
        int nO;
        nSi = nodes.n / 2.5;
        nO = nodes.n - nSi;
        xyzFile << nodes.n << std::endl;
        xyzFile << "# Element X Y Z  " << std::endl;
        xyzFile << std::fixed << std::showpoint << std::setprecision(6);
        for (int i = 0; i < nSi; ++i) {
            xyzFile << std::setw(10) << std::left << "Si" + std::to_string(i);
            xyzFile << std::setw(20) << std::left << nodes[i].crd[0];
            xyzFile << std::setw(20) << std::left << nodes[i].crd[1];
            xyzFile << std::setw(20) << std::left << "0.0" << std::endl;
        }
        for (int i = nSi; i < nodes.n; ++i) {
            xyzFile << std::setw(10) << std::left << "O" + std::to_string(i);
            xyzFile << std::setw(20) << std::left << nodes[i].crd[0];
            xyzFile << std::setw(20) << std::left << nodes[i].crd[1];
            xyzFile << std::setw(20) << std::left << "0.0" << std::endl;
        }
    }
    xyzFile.close();
}

void Network::writeBN(const std::string &prefix) {
    std::string newoutputPrefixfolder;
    std::string newoutputPrefixfile;

    newoutputPrefixfolder = prefix;
    newoutputPrefixfile = prefix;

    double BN_distance = 1.420 / 0.529177210903;
    double dim_x;
    double dim_y;

    dim_y = dimensions[1] * BN_distance;
    dim_x = dimensions[0] * BN_distance;
    std::ofstream dimFile(prefix + "_dimensions.dat",
                          std::ios::in | std::ios::trunc);

    dimFile << std::fixed << std::showpoint << std::setprecision(6);

    dimFile << dim_x << std::endl;
    dimFile << dim_y << std::endl;

    VecF<double> B(3);
    VecF<double> N(3);

    std::ofstream crysFile(prefix + "_BN_crys.crds",
                           std::ios::in | std::ios::trunc);
    crysFile << "T\nF\nF\nF\nF" << std::endl;
    crysFile << std::fixed << std::showpoint << std::setprecision(6);
    for (int i = 0; i < nodes.n; i = i + 2) {
        B[0] = nodes[i].crd[0] * BN_distance;
        B[1] = nodes[i].crd[1] * BN_distance;
        B[2] = 50.0;

        crysFile << std::setw(20) << std::left << B[0];
        crysFile << std::setw(20) << std::left << B[1];
        crysFile << std::setw(20) << std::left << B[2] << std::endl;
    }
    for (int i = 1; i < nodes.n; i = i + 2) {
        N[0] = nodes[i].crd[0] * BN_distance;
        N[1] = nodes[i].crd[1] * BN_distance;
        N[2] = 50.0;

        crysFile << std::setw(20) << std::left << N[0];
        crysFile << std::setw(20) << std::left << N[1];
        crysFile << std::setw(20) << std::left << N[2] << std::endl;
    }
    crysFile << std::setw(20) << std::left << 1.0;
    crysFile << std::setw(20) << std::left << 0.0;
    crysFile << std::setw(20) << std::left << 0.0 << std::endl;

    crysFile << std::setw(20) << std::left << 0.0;
    crysFile << std::setw(20) << std::left << 1.0;
    crysFile << std::setw(20) << std::left << 0.0 << std::endl;

    crysFile << std::setw(20) << std::left << 0.0;
    crysFile << std::setw(20) << std::left << 0.0;
    crysFile << std::setw(20) << std::left << 1.0 << std::endl;

    crysFile << std::setw(20) << std::left << dim_x << std::endl;
    crysFile << std::setw(20) << std::left << dim_y << std::endl;
    crysFile << std::setw(20) << std::left << 100.0 << std::endl;

    std::ofstream BNxyzFile(prefix + "_BN_crds.xyz",
                            std::ios::in | std::ios::trunc);

    BNxyzFile << nodes.n << std::endl;
    BNxyzFile << "# Element X Y Z  " << std::endl;
    BNxyzFile << std::fixed << std::showpoint << std::setprecision(6);
    for (int i = 0; i < nodes.n; i = i + 2) {
        B[0] = nodes[i].crd[0] * BN_distance;
        B[1] = nodes[i].crd[1] * BN_distance;
        B[2] = 50.0;
        BNxyzFile << std::setw(10) << std::left << "B" + std::to_string(int(i / 2));

        BNxyzFile << std::setw(20) << std::left << B[0];
        BNxyzFile << std::setw(20) << std::left << B[1];
        BNxyzFile << std::setw(20) << std::left << B[2] << std::endl;
    }
    for (int i = 1; i < nodes.n; i = i + 2) {
        N[0] = nodes[i].crd[0] * BN_distance;
        N[1] = nodes[i].crd[1] * BN_distance;
        N[2] = 50.0;
        BNxyzFile << std::setw(10) << std::left
                  << "N" + std::to_string(int((i - 1) / 2));

        BNxyzFile << std::setw(20) << std::left << N[0];
        BNxyzFile << std::setw(20) << std::left << N[1];
        BNxyzFile << std::setw(20) << std::left << N[2] << std::endl;
    }

    crysFile.close();
    BNxyzFile.close();
}
// Duplicated from procrystalline code

bool Network::r_ij(int i, int j, float cutoff) {

    float si_si_distance = 1.609 * sqrt((32.0 / 9.0));

    VecR<float> r_v(2, 2);
    VecR<float> r_v_original(2, 2);

    r_v[0] = (nodes[j].crd[0] - nodes[i].crd[0]) * si_si_distance;
    r_v[1] = (nodes[j].crd[1] - nodes[i].crd[1]) * si_si_distance;

    r_v_original[0] = r_v[0];
    r_v_original[1] = r_v[1];

    if (pow(r_v[0], 2) + pow(r_v[1], 2) > cutoff * cutoff) {
        for (int i = 0; i < 3; ++i) {
            r_v[0] = r_v_original[0] + (i - 1) * dimensions[0] * si_si_distance;
            for (int j = 0; j < 3; ++j) {
                r_v[1] = r_v_original[1] + (j - 1) * dimensions[1] * si_si_distance;
                if (pow(r_v[0], 2) + pow(r_v[1], 2) < cutoff * cutoff) {
                    return true;
                }
            }
        }
        return false;
    } else {
        return true;
    }
}

void Network::writeBilayerA(const std::string &prefix, float lj_cutoff) {
    std::cout << " ##### Starting Monolayer" << std::endl;
    float si_si_distance = 1.609 * sqrt((32.0 / 9.0));
    std::cout << "sisi    " << si_si_distance << std::endl;
    float o_o_distance = 1.609 * sqrt((8.0 / 3.0));
    float h = sin((19.5 / 180) * M_PI) * 1.609;

    VecF<double> diff;

    std::string newoutputPrefixfolder;
    std::string newoutputPrefixfile;

    newoutputPrefixfolder = prefix;
    newoutputPrefixfile = prefix;

    double dim_x;
    double dim_y;
    dim_y = dimensions[1] * si_si_distance;
    dim_x = dimensions[0] * si_si_distance;
    std::ofstream dimFile(newoutputPrefixfolder + "/dimensions.dat",
                          std::ios::in | std::ios::trunc);
    dimFile << std::fixed << std::showpoint << std::setprecision(6);

    dimFile << dim_x << std::endl;
    dimFile << dim_y << std::endl;

    std::ofstream crysFile(newoutputPrefixfolder + "/crys.crds",
                           std::ios::in | std::ios::trunc);
    std::ofstream harmpairsFile(newoutputPrefixfolder + "/harmpairs.dat",
                                std::ios::in | std::ios::trunc);
    std::ofstream bilayerxyzFile(newoutputPrefixfolder + "/bilayer.xyz",
                                 std::ios::in | std::ios::trunc);

    bilayerxyzFile << 6 * nodes.n << std::endl;
    bilayerxyzFile << "# Element X Y Z  " << std::endl;
    bilayerxyzFile << std::fixed << std::showpoint << std::setprecision(6);

    VecF<double> Si(3);
    VecF<double> O(3);
    VecF<int> Pair(2);

    int Si_O_harmpairs[2 * nodes.n][5];
    for (int i = 0; i < 2 * nodes.n; ++i) {
        for (int j = 0; j < 5; ++j) {
            Si_O_harmpairs[i][j] = 0;
        }
    }

    harmpairsFile << 10 * nodes.n * 2 << std::endl;
    crysFile << std::fixed << std::showpoint << std::setprecision(6);
    harmpairsFile << std::fixed << std::showpoint << std::setprecision(6);

    int atom_count = 1;
    std::cout << "##### Si crys" << std::endl;
    for (int i = 0; i < nodes.n; ++i) {
        Si_O_harmpairs[2 * i][0] = 2 * i + 1;
        Si_O_harmpairs[2 * i + 1][0] = 2 * i + 2;

        Si[0] = nodes[i].crd[0] * si_si_distance;
        Si[1] = nodes[i].crd[1] * si_si_distance;
        Si[2] = 5.0;

        crysFile << std::setw(20) << std::left << Si[0];
        crysFile << std::setw(20) << std::left << Si[1];
        crysFile << std::setw(20) << std::left << Si[2] << std::endl;

        bilayerxyzFile << std::setw(10) << std::left
                       << "Si" + std::to_string(2 * i);
        bilayerxyzFile << std::setw(20) << std::left << Si[0];
        bilayerxyzFile << std::setw(20) << std::left << Si[1];
        bilayerxyzFile << std::setw(20) << std::left << Si[2] << std::endl;

        atom_count += 1;
        Si[2] = 5.0 + 2.0 * 1.609;
        crysFile << std::setw(20) << std::left << Si[0];
        crysFile << std::setw(20) << std::left << Si[1];
        crysFile << std::setw(20) << std::left << Si[2] << std::endl;

        bilayerxyzFile << std::setw(10) << std::left
                       << "Si" + std::to_string(2 * i + 1);
        bilayerxyzFile << std::setw(20) << std::left << Si[0];
        bilayerxyzFile << std::setw(20) << std::left << Si[1];
        bilayerxyzFile << std::setw(20) << std::left << Si[2] << std::endl;

        atom_count += 1;
    }
    std::cout << "##### O crys" << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>> Axial " << std::endl;
    for (int i = 0; i < nodes.n; ++i) {
        O[0] = nodes[i].crd[0] * si_si_distance; // needs pbc work xoxo

        if (O[0] < 0)
            O[0] += dim_x;
        else if (O[0] > dim_x)
            O[0] -= dim_x;

        O[1] = nodes[i].crd[1] * si_si_distance;

        if (O[1] < 0)
            O[1] += dim_y;
        else if (O[1] > dim_y)
            O[1] -= dim_y;

        O[2] = 5 + 1.609;

        crysFile << std::setw(20) << std::left << O[0];
        crysFile << std::setw(20) << std::left << O[1];
        crysFile << std::setw(20) << std::left << O[2] << std::endl;
        bilayerxyzFile << std::setw(10) << std::left << "O" + std::to_string(i);
        bilayerxyzFile << std::setw(20) << std::left << O[0];
        bilayerxyzFile << std::setw(20) << std::left << O[1];
        bilayerxyzFile << std::setw(20) << std::left << O[2] << std::endl;

        Pair[0] = atom_count;
        Pair[1] = 2 * i + 1;

        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

        Si_O_harmpairs[2 * i][1] = atom_count;

        Pair[1] = 2 * i + 2;

        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

        Si_O_harmpairs[2 * i + 1][1] = atom_count;

        atom_count += 1;
    }

    for (int i = 0; i < nodes.n; ++i) {
        for (int j = 1; j < 5; j++) {
            nodes[i].oxyCnxs[j - 1] = Si_O_harmpairs[2 * i][j];
        }
        for (int j = 1; j < 5; j++) {
            nodes[i].oxyCnxs[j + 3] = Si_O_harmpairs[2 * i + 1][j];
        }
    }
    std::cout << ">>>>>>>>>> Equitorial" << std::endl;
    for (int i = 0; i < nodes.n; ++i) {
        for (int j = 0; j < nodes[i].netCnxs.n; ++j) {
            if (nodes[i].netCnxs[j] > i) {

                diff = (nodes[nodes[i].netCnxs[j]].crd - nodes[i].crd);
                diff[0] *= si_si_distance;
                diff[1] *= si_si_distance;
                float micx = dim_x / 2.0;
                float micy = dim_y / 2.0;

                if (diff[0] > micx)
                    diff[0] -= dim_x;
                else if (diff[0] < -micx)
                    diff[0] += dim_x;
                if (diff[1] > micy)
                    diff[1] -= dim_y;
                else if (diff[1] < -micy)
                    diff[1] += dim_y;

                O[0] = si_si_distance * (nodes[i].crd[0]) + diff[0] / 2.0;

                if (O[0] < 0)
                    O[0] += dim_x;
                else if (O[0] > dim_x)
                    O[0] -= dim_x;

                O[1] = si_si_distance * (nodes[i].crd[1]) + diff[1] / 2.;

                if (O[1] < 0)
                    O[1] += dim_y;
                else if (O[1] > dim_y)
                    O[1] -= dim_y;

                O[2] = 5 + 2 * 1.609 + h;
                // coordinate write
                crysFile << std::setw(20) << std::left << O[0];
                crysFile << std::setw(20) << std::left << O[1];
                crysFile << std::setw(20) << std::left << O[2] << std::endl;
                // bilayer write
                bilayerxyzFile << std::setw(10) << std::left
                               << "O" + std::to_string(nodes.n + 2 * i);
                bilayerxyzFile << std::setw(20) << std::left << O[0];
                bilayerxyzFile << std::setw(20) << std::left << O[1];
                bilayerxyzFile << std::setw(20) << std::left << O[2] << std::endl;

                Pair[0] = atom_count;
                Pair[1] = 2 * i + 1;

                harmpairsFile << std::setw(20) << std::left << Pair[0];
                harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

                Pair[1] = 2 * nodes[i].netCnxs[j] + 1;

                harmpairsFile << std::setw(20) << std::left << Pair[0];
                harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

                int k = 1;
                while (k < 5) {
                    if (Si_O_harmpairs[2 * i][k] == 0) {
                        Si_O_harmpairs[2 * i][k] = atom_count;
                        k += 5;
                    } else
                        ++k;
                }

                k = 1;
                while (k < 5) {
                    if (Si_O_harmpairs[2 * nodes[i].netCnxs[j]][k] == 0) {
                        Si_O_harmpairs[2 * nodes[i].netCnxs[j]][k] = atom_count;
                        k += 5;
                    } else
                        ++k;
                }
                atom_count += 1;

                O[2] = 5 - h;
                crysFile << std::setw(20) << std::left << O[0];
                crysFile << std::setw(20) << std::left << O[1];
                crysFile << std::setw(20) << std::left << O[2] << std::endl;

                bilayerxyzFile << std::setw(10) << std::left
                               << "O" + std::to_string(nodes.n + 2 * i + 1);
                bilayerxyzFile << std::setw(20) << std::left << O[0];
                bilayerxyzFile << std::setw(20) << std::left << O[1];
                bilayerxyzFile << std::setw(20) << std::left << O[2] << std::endl;

                Pair[0] = atom_count;
                Pair[1] = 2 * i + 2;
                harmpairsFile << std::setw(20) << std::left << Pair[0];
                harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

                Pair[1] = 2 * nodes[i].netCnxs[j] + 2;
                harmpairsFile << std::setw(20) << std::left << Pair[0];
                harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

                k = 1;
                while (k < 5) {
                    if (Si_O_harmpairs[2 * i + 1][k] == 0) {
                        Si_O_harmpairs[2 * i + 1][k] = atom_count;
                        k += 5;
                    } else
                        ++k;
                }
                k = 1;
                while (k < 5) {
                    if (Si_O_harmpairs[2 * nodes[i].netCnxs[j] + 1][k] == 0) {
                        Si_O_harmpairs[2 * nodes[i].netCnxs[j] + 1][k] = atom_count;
                        k += 5;
                    } else
                        ++k;
                }
                atom_count += 1;
            }
        }
    }

    std::cout << "##### O-O pairs" << std::endl;
    for (int i = 0; i < 2 * nodes.n; ++i) {
        Pair[0] = Si_O_harmpairs[i][1];
        Pair[1] = Si_O_harmpairs[i][2];
        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

        Pair[0] = Si_O_harmpairs[i][1];
        Pair[1] = Si_O_harmpairs[i][3];
        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

        Pair[0] = Si_O_harmpairs[i][1];
        Pair[1] = Si_O_harmpairs[i][4];
        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

        Pair[0] = Si_O_harmpairs[i][2];
        Pair[1] = Si_O_harmpairs[i][3];
        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

        Pair[0] = Si_O_harmpairs[i][2];
        Pair[1] = Si_O_harmpairs[i][4];
        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;

        Pair[0] = Si_O_harmpairs[i][3];
        Pair[1] = Si_O_harmpairs[i][4];
        harmpairsFile << std::setw(20) << std::left << Pair[0];
        harmpairsFile << std::setw(20) << std::left << Pair[1] << std::endl;
    }

    std::cout << "##### Tetrahedron" << std::endl;
    std::ofstream tetrahedraFile(newoutputPrefixfolder + "/tetrahedra.dat",
                                 std::ios::in | std::ios::trunc);
    VecF<int> Tetrahedron(5);
    for (int i = 0; i < 2 * nodes.n; ++i) {
        for (int j = 0; j < 5; j++) {
            Tetrahedron[j] = Si_O_harmpairs[i][j];
        }
        tetrahedraFile << std::setw(20) << std::left << Tetrahedron[0];
        tetrahedraFile << std::setw(20) << std::left << Tetrahedron[1];
        tetrahedraFile << std::setw(20) << std::left << Tetrahedron[2];
        tetrahedraFile << std::setw(20) << std::left << Tetrahedron[3];
        tetrahedraFile << std::setw(20) << std::left << Tetrahedron[4] << std::endl;
    }

    std::cout << "##### OPTIMISE_SILICA_INPUT" << std::endl;
    std::ofstream optimise_silica_File(newoutputPrefixfolder +
                                           "/optimise_silica.inpt",
                                       std::ios::in | std::ios::trunc);
    optimise_silica_File << "I/O" << std::endl;
    optimise_silica_File << "./crys.crds              input coordinates"
                         << std::endl;
    optimise_silica_File << "./harmpairs.dat              input harmonic pairs"
                         << std::endl;
    optimise_silica_File << "./lj_pairs.dat              input repulsive pairs"
                         << std::endl;
    optimise_silica_File << "./fixedz.dat              input fixed z atoms"
                         << std::endl;
    optimise_silica_File << "./test                              output prefix"
                         << std::endl;
    optimise_silica_File
        << "-----------------------------------------------------------"
        << std::endl;
    optimise_silica_File << "Restart Options" << std::endl;
    optimise_silica_File << "0               print restart file" << std::endl;
    optimise_silica_File << "0               restart run?" << std::endl;
    optimise_silica_File
        << "-----------------------------------------------------------"
        << std::endl;
    optimise_silica_File << "Sample Information" << std::endl;
    optimise_silica_File << int(nodes.n * 2) << std::endl;
    optimise_silica_File << "1               stretch x" << std::endl;
    optimise_silica_File << "1               stretch y" << std::endl;
    optimise_silica_File << "0               stretch z" << std::endl;
    optimise_silica_File << "45              central angle" << std::endl;
    optimise_silica_File << "0              scan angle" << std::endl;
    optimise_silica_File << double(dim_x) << std::endl;
    optimise_silica_File << double(dim_y) << std::endl;
    optimise_silica_File << "20              unit cell z" << std::endl;
    optimise_silica_File
        << "-----------------------------------------------------------"
        << std::endl;
    optimise_silica_File << "Geometry Optimisation" << std::endl;
    optimise_silica_File << "1                   resize(1/0)" << std::endl;
    optimise_silica_File
        << "1.300               starting volume (relative to reference)"
        << std::endl;
    optimise_silica_File
        << "0.900               final volume (relative to reference)"
        << std::endl;
    optimise_silica_File << "5                   number of volumes to analyse"
                         << std::endl;
    optimise_silica_File << "1                  samples per volume" << std::endl;
    optimise_silica_File << double(dim_x) << std::endl;
    optimise_silica_File << double(dim_y) << std::endl;
    optimise_silica_File << "20                  unit cell z reference"
                         << std::endl;
    optimise_silica_File
        << "10000000             max steps iterations steepest descent (per area)"
        << std::endl;
    optimise_silica_File << "0.5                 Armijo backtracking constant"
                         << std::endl;
    optimise_silica_File << "1e-9                convergence tolerance"
                         << std::endl;
    optimise_silica_File
        << "-----------------------------------------------------------"
        << std::endl;
    optimise_silica_File << "Potential Model" << std::endl;
    optimise_silica_File
        << "1                   turn on harmonic interactions (1/0)" << std::endl;
    optimise_silica_File << "1.609               harmonic Si-O r0" << std::endl;
    optimise_silica_File << "1                   harmonic Si-O k" << std::endl;
    optimise_silica_File << "1                   harmonic O-O k" << std::endl;
    optimise_silica_File << "0                   turn on SiSi harmonics (1/0)"
                         << std::endl;
    optimise_silica_File << "0                   harmonic Si-Si r0" << std::endl;
    optimise_silica_File << "0                   harmonic Si-Si k" << std::endl;
    optimise_silica_File << "1                   turn on 24-12 repulsions (1/0)"
                         << std::endl;
    optimise_silica_File << "1.7                 repulsive r0" << std::endl;
    optimise_silica_File << "0.25                repulsive k" << std::endl;
    optimise_silica_File << "0                   turn on z fixing (1/0)"
                         << std::endl;
    optimise_silica_File << "1                   turn on multicore optimisation"
                         << std::endl;
    optimise_silica_File << "0                   Parallelise Sample" << std::endl;
    optimise_silica_File << "1                   Parallelise Area" << std::endl;
    optimise_silica_File << "1                   number of cores" << std::endl;
    optimise_silica_File << "0                   turn on cuda" << std::endl;
    optimise_silica_File
        << "-----------------------------------------------------------"
        << std::endl;

    crysFile.close();
    harmpairsFile.close();
    bilayerxyzFile.close();

    tetrahedraFile.close();

    optimise_silica_File.close();

    std::cout << "Beginning Lj File" << std::endl;
    std::ofstream ljFile(newoutputPrefixfolder + "/lj_pairs.dat",
                         std::ios::in | std::ios::trunc);

    VecR<int> lj_i(100000, 100000);
    VecR<int> lj_j(100000, 100000);

    std::cout << "LJ CUTOFF : " << lj_cutoff << std::endl;

    int lj_count = 0;

    for (int i = 0; i < nodes.n; ++i) {
        lj_i[lj_count] = (2 * i + 1);
        lj_j[lj_count] = (2 * i + 2);
        lj_count += 1;
    }

    for (int i = 0; i < nodes.n - 1; ++i) {
        for (int j = i + 1; j < nodes.n; ++j) {
            if (r_ij(i, j, lj_cutoff) == true) {

                lj_i[lj_count] = (2 * i + 1);
                lj_j[lj_count] = (2 * j + 1);
                lj_count += 1;

                lj_i[lj_count] = (2 * i + 2);
                lj_j[lj_count] = (2 * j + 2);
                lj_count += 1;

                lj_i[lj_count] = (2 * i + 1);
                lj_j[lj_count] = (2 * j + 2);
                lj_count += 1;

                lj_i[lj_count] = (2 * i + 2);
                lj_j[lj_count] = (2 * j + 1);
                lj_count += 1;
            }
        }
    }

    ljFile << lj_count << std::endl;
    ljFile << std::fixed << std::showpoint << std::setprecision(6);

    for (int i = 0; i < lj_count; ++i) {
        ljFile << std::setw(20) << std::left << lj_i[i];
        ljFile << std::setw(20) << std::left << lj_j[i] << std::endl;
    }

    ljFile.close();
}

// Create a Triangle raft
// The idea here is that these will more accurately reflect silica potentials
// This is not quite as easy as originally thought, as it requires additional
// layers to the current code As such, the triangle raft will run ALONGSIDE the
// graphitic like potential, with bond switches occuring in both The main
// difference is that energy will be calculated for this lattice, rather than
// networkA However, all changes must be mirrored in all networks

Network::Network(VecR<Node> nodesA, VecF<double> pbA, VecF<double> rpbA,
                 int maxNetCnxsA, int maxDualCnxsA)
    : dimensions(pbA), rpb(rpbA), maxNetCnxs(maxNetCnxsA * 2),
      maxDualCnxs(maxDualCnxsA) {
    // Taking Nodes and Periodic details from networkA to build new network
    int nNodes;
    int OxygenAtom = nodesA.n;
    int atom0;
    int atom1;
    VecF<double> crds(2);
    VecF<double> crds1(2);
    VecF<double> crds2(2);
    VecF<double> v(2);

    nNodes = nodesA.n * 5 / 2;
    geometryCode = "2DEtr";

    // Get Graphene -> Triangle raft conversion array;
    // Convert Graphene Metal atoms -> Triangle raft atoms
    // Add oxygens
    nodes = VecR<Node>(0, nNodes);
    for (int i = 0; i < nodesA.n; ++i) {
        Node node(i, maxNetCnxs, maxDualCnxs, 0);
        nodes.addValue(node);
        for (int j = 0; j < nodesA[i].netCnxs.n; ++j) {
            nodes[i].netCnxs.addValue(nodesA[i].netCnxs[j]);
        }
        for (int j = 0; j < nodesA[i].dualCnxs.n; ++j) {
            nodes[i].dualCnxs.addValue(nodesA[i].dualCnxs[j]);
        }
        nodes[i].crd = nodesA[i].crd * SiScaling;
    }
    for (int i = 0; i < nodesA.n; ++i) {
        atom0 = i;
        crds1[0] = nodes[atom0].crd[0];
        crds1[1] = nodes[atom0].crd[1];

        for (int j = 0; j < nodes[i].netCnxs.n; ++j) {
            atom1 = nodes[i].netCnxs[j];
            crds2[0] = nodes[atom1].crd[0];
            crds2[1] = nodes[atom1].crd[1];

            if (atom1 > i && atom1 < nodesA.n) {

                v[0] = crds1[0] - crds2[0];
                if (v[0] > dimensions[0] / 2.0) {
                    v[0] -= dimensions[0];
                } else if (v[0] < -dimensions[0] / 2.0) {
                    v[0] += dimensions[0];
                }
                v[1] = crds1[1] - crds2[1];
                if (v[1] > dimensions[1] / 2.0) {
                    v[1] -= dimensions[1];
                } else if (v[1] < -dimensions[1] / 2.0) {
                    v[1] += dimensions[1];
                }

                // Fit oxygens between these two Si atoms.
                crds[0] = crds2[0] + v[0] / 2.0;
                crds[1] = crds2[1] + v[1] / 2.0;

                if (crds[0] > dimensions[0]) {
                    crds[0] -= dimensions[0];
                } else if (crds[0] < -dimensions[0]) {
                    crds[0] += dimensions[0];
                }
                if (crds[1] > dimensions[1]) {
                    crds[1] -= dimensions[1];
                } else if (crds[1] < -dimensions[1]) {
                    crds[1] += dimensions[1];
                }

                Node node(OxygenAtom, maxNetCnxs * 2, maxDualCnxs, 0);
                nodes.addValue(node);

                nodes[OxygenAtom].crd = crds;
                nodes[atom0].netCnxs.swapValue(atom1, OxygenAtom, true);
                nodes[atom1].netCnxs.swapValue(atom0, OxygenAtom, true);
                nodes[OxygenAtom].netCnxs.addValue(atom0);
                nodes[OxygenAtom].netCnxs.addValue(atom1);
                OxygenAtom += 1;
            }
        }
    }
    for (int i = 0; i < nodesA.n; ++i) {
        int atom1 = nodes[i].netCnxs[0];
        int atom2 = nodes[i].netCnxs[1];
        int atom3 = nodes[i].netCnxs[2];

        nodes[atom1].netCnxs.addValue(atom2);
        nodes[atom2].netCnxs.addValue(atom1);

        nodes[atom1].netCnxs.addValue(atom3);
        nodes[atom3].netCnxs.addValue(atom1);

        nodes[atom2].netCnxs.addValue(atom3);
        nodes[atom3].netCnxs.addValue(atom2);
    }
}

int Network::getMaxCnxs() {
    int maxNetCnxs = 0;
    for (int i = 0; i < nodes.n; ++i) {
        if (nodes[i].netCnxs.n > maxNetCnxs) {
            maxNetCnxs = nodes[i].netCnxs.n;
        }
    }
    return maxNetCnxs;
}

int Network::getMinCnxs() {
    int minNetCnxs = 100000;
    for (int i = 0; i < nodes.n; ++i) {
        if (nodes[i].netCnxs.n < minNetCnxs) {
            minNetCnxs = nodes[i].netCnxs.n;
        }
    }
    return minNetCnxs;
}