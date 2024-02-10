#include "linked_network.h"
#include <algorithm>
#include <chrono>
#include <ctime>
#include <iostream>
#include <omp.h>
#include <unistd.h>

// Default constructor
LinkedNetwork::LinkedNetwork() = default;

/**
 * @brief Construct a hexagonal linked network from scratch using netmc.inpt
 * parameters
 * @param numRings the number of nodes in lattice A
 * @param logger the logger object
 */
LinkedNetwork::LinkedNetwork(int numRings, LoggerPtr logger)
    : networkB(numRings), minBCnxs(6), maxBCnxs(6), minACnxs(3), maxACnxs(3), centreCoords(2) {

    networkA = networkB.constructDual(maxACnxs);
    networkA.maxNetCnxs = maxACnxs;
    rescale(sqrt(3));
    networkB.gen2ndOrderConnections(networkA);
    dimensions = networkA.dimensions;
    centreCoords[0] = dimensions[0] / 2;
    centreCoords[1] = dimensions[1] / 2;
}

/**
 * @brief Construct by loading networks from files
 * @param inputData the input data object
 * @param logger the logger object
 */
LinkedNetwork::LinkedNetwork(InputData &inputData, LoggerPtr logger) : centreCoords(2) {
    std::string prefix = inputData.inputFolder + '/' + inputData.inputFilePrefix;

    networkA = Network(prefix + "_A", inputData.maxRingSize, inputData.maxRingSize, logger);
    minACnxs = networkA.getMinCnxs();
    maxACnxs = networkA.getMaxCnxs();

    networkB = Network(prefix + "_B", inputData.maxRingSize, inputData.maxRingSize, logger);
    minBCnxs = networkB.getMinCnxs();
    maxBCnxs = networkB.getMaxCnxs();

    dimensions = networkA.dimensions;
    centreCoords[0] = dimensions[0] / 2;
    centreCoords[1] = dimensions[1] / 2;

    if (inputData.minCoordination <= minACnxs) {
        minACnxs = inputData.minCoordination;
    } else {
        logger->warn("Initial base network does not fit within allowed min node "
                     "coordination numbers, input file: {} vs NetMC files: {}",
                     inputData.minCoordination, minACnxs);
    }
    if (inputData.maxCoordination >= maxACnxs) {
        maxACnxs = inputData.maxCoordination;
    } else {
        logger->warn("Initial base network does not fit within allowed max node "
                     "coordination numbers, input file: {} vs NetMC files: {}",
                     inputData.maxCoordination, maxACnxs);
    }
    if (inputData.minRingSize <= minBCnxs) {
        minBCnxs = inputData.minRingSize;
    } else {
        logger->warn("Initial ring network does not fit within allowed min ring "
                     "coordination numbers, input file: {} vs NetMC files: {}",
                     inputData.minRingSize, minBCnxs);
    }
    if (inputData.maxRingSize >= maxBCnxs) {
        maxBCnxs = inputData.maxRingSize;
    } else {
        logger->warn("Initial ring network does not fit within allowed max ring "
                     "coordination numbers, input file: {} vs NetMC files: {}",
                     inputData.maxRingSize, maxBCnxs);
    }
    logger->info("Creating network T...");
    if (inputData.isRestartUsingLAMMPSObjectsEnabled && inputData.isTriangleRaftEnabled) {
        networkT = Network(prefix + "_Si2O3", inputData.maxRingSize, inputData.maxRingSize, logger);
    } else {
        networkT = Network(networkA.nodes, networkA.dimensions, networkA.rpb, maxACnxs, maxACnxs);
    }

    //  Linked Network can contain lammps objects for minimisation
    //  These will not necessarily have the same coordinate systems, bond orders,
    //  etc etc.
    if (inputData.isSimpleGrapheneEnabled)
        SimpleGraphene = LammpsObject("Si", inputData.inputFolder, logger);
    if (inputData.isTersoffGrapheneEnabled)
        TersoffGraphene = LammpsObject("C", inputData.inputFolder, logger);
    if (inputData.isTriangleRaftEnabled)
        Triangle_Raft = LammpsObject("Si2O3", inputData.inputFolder, logger);
    if (inputData.isBilayerEnabled)
        Bilayer = LammpsObject("SiO2", inputData.inputFolder, logger);
    if (inputData.isBNEnabled)
        BN = LammpsObject("BN", inputData.inputFolder, logger);

    if (inputData.isFixRingsEnabled) {
        findFixedRings(inputData.isFixRingsEnabled, prefix, logger);
        findFixedNodes();
    }
}

void LinkedNetwork::findFixedNodes() {
    fixedNodes.setSize(0);
    for (int i = 0; i < fixedRings.n; ++i) {
        for (int j = 0; j < networkB.nodes[fixedRings[i]].dualCnxs.n; ++j) {
            if (fixedNodes.contains(networkB.nodes[fixedRings[i]].dualCnxs[j]) == false)
                fixedNodes.addValue(networkB.nodes[fixedRings[i]].dualCnxs[j]);
        }
    }
}

/**
 * @brief read the fixed_rings.dat file and store the fixed ring IDs in a vector
 * @param isFixedRingsEnabled boolean to enable or disable fixed rings
 * @param filename the name of the input file
 * @param logger the logger object
 */
void LinkedNetwork::findFixedRings(bool isFixedRingsEnabled,
                                   std::string filename, LoggerPtr logger) {
    // Format of the fixed_rings.dat file has changed to not include the number of
    // fixed rings on the first line, but simply the IDs of each fixed ring line
    // by line.
    if (isFixedRingsEnabled) {
        // Open the file
        std::ifstream fixedRingsFile(filename + ".dat", std::ios::in);
        if (!fixedRingsFile.is_open()) {
            logger->warn("Failed to open file: {}.dat, setting number of fixed rings to 0.", filename);
            fixedRings.setSize(0);
            return;
        }
        std::string line;
        std::string ringList = "";
        int numFixedRings = 0;
        // Read the file line by line
        while (getline(fixedRingsFile, line)) {
            // Convert the line to an integer and store it in the vector
            std::istringstream(line) >> fixedRings[numFixedRings];
            ringList += std::to_string(fixedRings[numFixedRings]) + " ";
            numFixedRings++;
        }
        fixedRings.setSize(numFixedRings);
        logger->info("Number of fixed rings: {}", numFixedRings);
        logger->info("Fixed rings: {}", ringList);
    } else {
        logger->info("Fixed rings disabled, setting number of fixed rings to 0.");
        fixedRings.setSize(0);
    }
}

void LinkedNetwork::pushPrefix(std::string_view prefixInArg,
                               std::string_view prefixOutArg) {
    prefixIn = prefixInArg;
    prefixOut = prefixOutArg;
}

// Set up potential model with single angle and bond parameter set
void LinkedNetwork::initialisePotentialModel(double harmonicAngleForceConstant,
                                             double harmonicBondForceConstant,
                                             double harmonicGeometryConstraint,
                                             bool isMaintainConvexityEnabledArg,
                                             LoggerPtr logger) {
    // Make copy of lattice A coordinates
    logger->info("Copying lattice A coordinates");
    crds = VecF<double>(2 * networkA.nodes.n);
    for (int i = 0; i < networkA.nodes.n; ++i) {
        crds[2 * i] = networkA.nodes[i].crd[0];
        crds[2 * i + 1] = networkA.nodes[i].crd[1];
    }

    // Initialise potential model parameters
    // Angle parameters
    logger->info("Initialising potential model parameters");
    potParamsA = VecF<double>(6); // for 3 and 4 coordinate
    potParamsA[0] = harmonicAngleForceConstant;
    potParamsA[1] = cos(2 * M_PI / 3.0);
    potParamsA[2] = harmonicBondForceConstant;
    potParamsA[3] = cos(2 * M_PI / 4.0);

    // Bond parameters
    logger->info("Initialising bond parameters");
    potParamsB = VecF<double>(3);
    potParamsB[0] = harmonicBondForceConstant;
    potParamsB[1] = 1.0;

    // Geometry constraint parameters
    logger->info("Initialising geometry constraint parameters");
    potParamsC = VecF<double>(2);
    potParamsC[0] =
        harmonicGeometryConstraint; // k, r0 updated through optimal projection
    isMaintainConvexityEnabled = isMaintainConvexityEnabledArg;
}

// Set up geometry optimisation parameters
void LinkedNetwork::initialiseGeometryOpt(int iterations, double tau,
                                          double tolerance, int localExtent) {
    goptParamsA = VecF<int>(2);
    goptParamsA[0] = iterations;
    goptParamsA[1] = localExtent;
    goptParamsB = VecF<double>(2);
    goptParamsB[0] = tau;
    goptParamsB[1] = tolerance;
}

// Set up monte carlo and random number generators
void LinkedNetwork::initialiseMonteCarlo(const Network &network,
                                         double temperature, LoggerPtr logger,
                                         int seed) {
    logger->info("Collecting global energy...");
    double energy =
        globalPotentialEnergy(false, isMaintainConvexityEnabled, network, logger);
    logger->info("Global optimisation energy : {}", energy);
    mc = Metropolis(seed, temperature, energy);
    mtGen.seed(seed);
}

// Rescale lattice dimensions
void LinkedNetwork::rescale(double scaleFactor) {
    networkA.rescale(scaleFactor);
    networkB.rescale(scaleFactor);
}

// Constants for connection types
const int CNX_TYPE_33 = 33;
const int CNX_TYPE_44 = 44;
const int CNX_TYPE_43 = 43;
/**
 * @brief Assigns values to a, b, and cnxType based on the connections of the two random nodes.
 * @param a Reference to ID of the first node in the selected connection.
 * @param b Reference to ID of the second node in the selected connection.
 * @param cnxType Reference to the type of the selected connection. 33, 44, or 43.
 * @param randNode ID of the first random node.
 * @param randNodeConnection ID of the second random node.
 * @param randNodeCoordination Coordination of the first random node.
 * @param randNodeConnectionCoordination Coordination of the second random node.
 * @param logger The logger object.
 */
void LinkedNetwork::assignValues(int &a, int &b, int &cnxType, int randNode, int randNodeConnection,
                                 int randNodeCoordination, int randNodeConnectionCoordination) {
    if (randNodeCoordination == 3 && randNodeConnectionCoordination == 3) {
        cnxType = CNX_TYPE_33;
        a = randNode;
        b = randNodeConnection;
    } else if (randNodeCoordination == 4 && randNodeConnectionCoordination == 4) {
        cnxType = CNX_TYPE_44;
        a = randNode;
        b = randNodeConnection;
    } else if ((randNodeCoordination == 3 && randNodeConnectionCoordination == 4) ||
               (randNodeCoordination == 4 && randNodeConnectionCoordination == 3)) {
        cnxType = CNX_TYPE_43;
        a = (randNodeCoordination == 4) ? randNode : randNodeConnection;
        b = (randNodeCoordination == 4) ? randNodeConnection : randNode;
    } else {
        std::ostringstream oss;
        oss << "Nodes in random connection have unsupported coordinations: " << randNodeCoordination
            << " and " << randNodeConnectionCoordination;
        throw std::runtime_error(oss.str());
    }
}

/**
 * @brief Assigns values to baseNode1, baseNode2, ringNode1, ringNode2 by chosing a random connection between nodes in
 * lattice A with a weighted probability based on the distance of each node from the centre of the network.
 * @param baseNode1 Reference to ID of the first base node to be assigned.
 * @param baseNode2 Reference to ID of the second base node to be assigned.
 * @param ringNode1 Reference to ID of the first ring node to be assigned.
 * @param ringNode2 Reference to ID of the second ring node to be assigned.
 * @param gen Reference to a Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
 * @param selectionType The type of selection to be used. Exponential decay or random.
 * @param logger The logger object.
 * @return The type of the selected connection. 33, 34 or 44.
 * @throw std::runtime_error if the nodes in the random connection have coordinations other than 3 or 4.
 */
int LinkedNetwork::pickRandomConnection(int &baseNode1, int &baseNode2, int &ringNode1, int &ringNode2,
                                        std::mt19937 &gen, SelectionType selectionType, LoggerPtr logger) {

    /*
     * 3-3 coordination connection
     * a,b,c,d,e,f are nodes in lattice A
     * u,v,w,x are nodes in lattice B
     * E       F            r1
     *  \     /           / | \
     *   b1--b2          W  |  X
     *  /     \           \ | /
     * C       D            r2
     */

    std::vector<double> weights(networkA.nodes.n);
    if (selectionType == SelectionType::EXPONENTIAL_DECAY) {
        double scale = 1.0; // Adjust this value based on your needs
        for (int i = 0; i < networkA.nodes.n; ++i) {
            double distance = networkA.nodes[i].distanceFrom(centreCoords);
            weights[i] = std::exp(-distance * scale);
        }

        // Normalize weights
        double total = std::accumulate(weights.begin(), weights.end(), 0.0);
        for (double &weight : weights) {
            weight /= total;
        }
    } else { // SelectionType::RANDOM
        std::fill(weights.begin(), weights.end(), 1.0);
    }

    // Create a discrete distribution based on the weights
    std::discrete_distribution<> distribution(weights.begin(), weights.end());
    int cnxType;
    bool pickingAcceptableRing = true;
    while (pickingAcceptableRing) {
        int randNode = distribution(gen);
        int randNodeCoordination = networkA.nodes[randNode].netCnxs.n;
        std::uniform_int_distribution<int> randomCnx(0, randNodeCoordination - 1);
        int randNodeConnection = networkA.nodes[randNode].netCnxs[randomCnx(gen)];
        int randNodeConnectionCoordination = networkA.nodes[randNodeConnection].netCnxs.n;

        // Use helper function to assign a, b, and cnxType
        try {
            assignValues(baseNode1, baseNode2, cnxType, randNode, randNodeConnection, randNodeCoordination,
                         randNodeConnectionCoordination);
        } catch (std::runtime_error &e) {
            logger->critical(e.what());
            throw;
        }
        // Two connected nodes should always share two ring nodes.
        if (VecR<int> commonRings = vCommonValues(networkA.nodes[baseNode1].dualCnxs, networkA.nodes[baseNode2].dualCnxs);
            commonRings.n == 2) {
            // Randomly assign u and v to those two common ring nodes
            std::uniform_int_distribution<int> randomDirection(0, 1);
            int randIndex = randomDirection(gen);
            ringNode1 = commonRings[randIndex];
            ringNode2 = commonRings[1 - randIndex];
        } else {
            wrapCoordinates();
            syncCoordinates();
            write("debug");
            logger->critical("Selected random connection does not share two ring nodes");
            throw std::runtime_error("Selected random connection does not share two ring nodes");
        }
        // Check that the base nodes are not a member of a fixed ring.
        if (!fixedNodes.contains(baseNode1) && !fixedNodes.contains(baseNode2)) {
            // If so, break and return the connection type
            pickingAcceptableRing = false;
        }
    }
    return cnxType;
}

/**
 * @brief Generate all ids of nodes in lattices A and B needed for switch move,
 * only for 3/4 coordinate nodes
 * @param cnxType the type of connection to be generated
 * @param switchIdsA the ids of the nodes in lattice A
 * @param switchIdsB the ids of the nodes in lattice B
 * @param switchIdsT the ids of the nodes in lattice T
 * @param a the id of the first node in lattice A
 * @param b the id of the second node in lattice A
 * @param u the id of the first node in lattice B
 * @param v the id of the second node in lattice B
 * @return true if the switch move is possible, false otherwise
 */
bool LinkedNetwork::generateSwitchIDs(int cnxType, VecF<int> &switchIdsA, VecF<int> &switchIdsB, VecF<int> &switchIdsT,
                                      int a, int b, int u, int v, LoggerPtr logger) {

    // lots of error checking to remove any potential pathological cases
    if (a == b || u == v)
        return false;
    if (cnxType == 33) {
        /* Switch connectivities in lattice and dual
         * 3-3 coordination connection
         * a,b,c,d,e,f are nodes in lattice A
         * u,v,w,x are nodes in lattice B
         *  E      F            V
         *   \    /           / | \
         *   A---B           W  |  X
         *  /     \           \ | /
         * C       D            U
         */

        int errorFlag = 0;
        VecR<int> common;
        VecR<int> common1;

        int c = findCommonConnection(a, u, b);
        int d = findCommonConnection(b, u, a);
        int e = findCommonConnection(a, v, b);
        int f = findCommonConnection(b, v, a);
        int w = findAssociatedNodeAA(a, c, u);
        int x = findAssociatedNodeAA(b, d, u);
        int alpha = 0;
        int beta = 0;
        int gamma = 0;
        int delta = 0;
        int eta = 0;
        for (int i = 0; i < networkT.nodes[a].netCnxs.n; ++i) {
            for (int j = 0; j < networkT.nodes[b].netCnxs.n; ++j) {
                if (networkT.nodes[a].netCnxs[i] == networkT.nodes[b].netCnxs[j]) {
                    alpha = networkT.nodes[a].netCnxs[i];
                }
            }
            for (int j = 0; j < networkT.nodes[c].netCnxs.n; ++j) {
                if (networkT.nodes[a].netCnxs[i] == networkT.nodes[c].netCnxs[j]) {
                    beta = networkT.nodes[a].netCnxs[i];
                }
            }
            for (int j = 0; j < networkT.nodes[e].netCnxs.n; ++j) {
                if (networkT.nodes[a].netCnxs[i] == networkT.nodes[e].netCnxs[j]) {
                    gamma = networkT.nodes[a].netCnxs[i];
                }
            }
        }
        for (int i = 0; i < networkT.nodes[b].netCnxs.n; ++i) {
            for (int j = 0; j < networkT.nodes[d].netCnxs.n; ++j) {
                if (networkT.nodes[b].netCnxs[i] == networkT.nodes[d].netCnxs[j]) {
                    delta = networkT.nodes[b].netCnxs[i];
                }
            }
            for (int j = 0; j < networkT.nodes[f].netCnxs.n; ++j) {
                if (networkT.nodes[b].netCnxs[i] == networkT.nodes[f].netCnxs[j]) {
                    eta = networkT.nodes[b].netCnxs[i];
                }
            }
        }
        if (alpha == 0) {
            logger->error("alpha broken");
        }
        if (beta == 0) {
            logger->error("beta broken");
        }
        if (gamma == 0) {
            logger->error("gamma broken");
        }
        if (delta == 0) {
            logger->error("delta broken");
        }
        if (eta == 0) {
            logger->error("eta broken");
        }

        // Additional error checking
        if (c == d || e == f)
            errorFlag = 6; // can simply be triangle edge sharing pair (not an error)
        // Prevent rings having only two or fewer neighbours
        VecR<int> vCnxs = vUnique(networkB.nodes[v].netCnxs);
        if (vCnxs.n <= 3)
            errorFlag = 10;
        VecR<int> uCnxs = vUnique(networkB.nodes[u].netCnxs);
        uCnxs.delValue(v);
        if (uCnxs.n <= 2) {
            for (int i = 0; i < uCnxs.n; ++i) {
                if (uCnxs[i] == x) {
                    errorFlag = 10;
                    break;
                }
            }
        }
        if (errorFlag != 0)
            return false;
        // check move will not violate dual connectivity limits
        if (networkB.nodes[u].netCnxs.n == minBCnxs ||
            networkB.nodes[v].netCnxs.n == minBCnxs ||
            networkB.nodes[w].netCnxs.n == maxBCnxs ||
            networkB.nodes[x].netCnxs.n == maxBCnxs)
            return false;

        switchIdsA = VecF<int>(6);
        switchIdsA[0] = a;
        switchIdsA[1] = b;
        switchIdsA[2] = c;
        switchIdsA[3] = d;
        switchIdsA[4] = e;
        switchIdsA[5] = f;

        switchIdsB = VecF<int>(4);
        switchIdsB[0] = u;
        switchIdsB[1] = v;
        switchIdsB[2] = w;
        switchIdsB[3] = x;

        switchIdsT = VecF<int>(11);
        switchIdsT[0] = a;
        switchIdsT[1] = b;
        switchIdsT[2] = c;
        switchIdsT[3] = d;
        switchIdsT[4] = e;
        switchIdsT[5] = f;
        switchIdsT[6] = alpha;
        switchIdsT[7] = beta;
        switchIdsT[8] = gamma;
        switchIdsT[9] = delta;
        switchIdsT[10] = eta;
        return true;
    } else if (cnxType == 44) {
        /* Switch connectivities in lattice and dual
         * 4-4 coordination connection
         * a,b,c,d,e,f,g,h are nodes in lattice A
         * u,v,w,x,y,z are nodes in lattice B
         * a-b, a-c, a-e, a-g
         * b-a, b-d, b-f, b-h
         * a-b share u-v
         * c-a-b-d share u
         * e-a-b-f share v
         * g-a-e share w
         * d-b-h share x
         * g-a-c share y
         * f-b-h share z
         * u-v, u-y, u-x
         * v-u, v-w, v-z
         * w-y, x-z*/
        int errorFlag = 0;
        int c;
        int d;
        int e;
        int f;
        int g;
        int h;
        int w;
        int x;
        int y;
        int z;

        VecR<int> common;
        VecR<int> common1;
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(b);
        if (common.n != 1)
            errorFlag = 1;
        c = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(a);
        if (common.n != 1)
            errorFlag = 2;
        d = common[0];
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(b);
        if (common.n != 1)
            errorFlag = 3;
        e = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(a);
        if (common.n != 1)
            errorFlag = 4;
        f = common[0];

        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[e].dualCnxs);
        common.delValue(v);
        if (common.n != 1)
            errorFlag = 5;
        w = common[0];
        common = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[d].dualCnxs);
        common.delValue(u);
        if (common.n != 1)
            errorFlag = 6;
        x = common[0];
        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[c].dualCnxs);
        common.delValue(u);
        if (common.n != 1)
            errorFlag = 7;
        y = common[0];
        common = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[f].dualCnxs);
        common.delValue(v);
        if (common.n != 1)
            errorFlag = 8;
        z = common[0];

        common = networkA.nodes[a].netCnxs;
        common.delValue(b);
        common.delValue(c);
        common.delValue(e);
        if (common.n != 1)
            errorFlag = 9;
        g = common[0];
        common = networkA.nodes[b].netCnxs;
        common.delValue(a);
        common.delValue(d);
        common.delValue(f);
        if (common.n != 1)
            errorFlag = 10;
        h = common[0];

        // Additional error checking including preventing two nodes connecting
        // multiple times
        if (c == d || e == f)
            errorFlag = 6;
        if (vContains(networkB.nodes[w].netCnxs, x))
            errorFlag = 7;
        if (vContains(networkB.nodes[y].netCnxs, x))
            errorFlag = 7;
        if (vContains(networkB.nodes[y].auxCnxs, x))
            errorFlag = 7;
        if (vContains(networkB.nodes[z].netCnxs, w))
            errorFlag = 7;
        if (vContains(networkB.nodes[z].auxCnxs, w))
            errorFlag = 7;

        if (errorFlag != 0)
            return false;
        // check move will not violate dual connectivity limits
        if (networkB.nodes[u].netCnxs.n == minBCnxs ||
            networkB.nodes[v].netCnxs.n == minBCnxs ||
            networkB.nodes[w].netCnxs.n == maxBCnxs ||
            networkB.nodes[x].netCnxs.n == maxBCnxs)
            return false;

        switchIdsA = VecF<int>(8);
        switchIdsA[0] = a;
        switchIdsA[1] = b;
        switchIdsA[2] = c;
        switchIdsA[3] = d;
        switchIdsA[4] = e;
        switchIdsA[5] = f;
        switchIdsA[6] = g;
        switchIdsA[7] = h;

        switchIdsB = VecF<int>(6);
        switchIdsB[0] = u;
        switchIdsB[1] = v;
        switchIdsB[2] = w;
        switchIdsB[3] = x;
        switchIdsB[4] = y;
        switchIdsB[5] = z;
        return true;
    } else if (cnxType == 43) {
        /* Switch connectivities in lattice and dual
         * 4-3 coordination connection
         * a,b,c,d,e,f,g are nodes in lattice A
         * u,v,w,x,y are nodes in lattice B
         * a-b, a-c, a-e, a-g
         * b-a, b-d, b-f
         * a-b share u-v
         * c-a-b-d share u
         * e-a-b-f share v
         * g-a-e share w
         * d-b-f share x
         * g-a-c share y
         * u-v, u-y, u-x
         * v-u, v-w, v-x
         * w-y */
        int errorFlag = 0;
        int c;
        int d;
        int e;
        int f;
        int g;
        int w;
        int x;
        int y;

        VecR<int> common;
        VecR<int> common1;
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(b);
        if (common.n != 1)
            errorFlag = 1;
        c = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[u].dualCnxs);
        common.delValue(a);
        if (common.n != 1)
            errorFlag = 2;
        d = common[0];
        common = vCommonValues(networkA.nodes[a].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(b);
        if (common.n != 1)
            errorFlag = 3;
        e = common[0];
        common = vCommonValues(networkA.nodes[b].netCnxs, networkB.nodes[v].dualCnxs);
        common.delValue(a);
        if (common.n != 1)
            errorFlag = 4;
        f = common[0];

        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[e].dualCnxs);
        common.delValue(v);
        if (common.n != 1)
            errorFlag = 5;
        w = common[0];
        common = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[d].dualCnxs);
        common1 = vCommonValues(networkA.nodes[b].dualCnxs, networkA.nodes[f].dualCnxs);
        common.delValue(u);
        common1.delValue(v);
        if (common.n != 1 || common1.n != 1 || common[0] != common1[0])
            errorFlag = 5;
        x = common[0];
        common = vCommonValues(networkA.nodes[a].dualCnxs, networkA.nodes[c].dualCnxs);
        common.delValue(u);
        if (common.n != 1)
            errorFlag = 7;
        y = common[0];

        common = networkA.nodes[a].netCnxs;
        common.delValue(b);
        common.delValue(c);
        common.delValue(e);
        if (common.n != 1)
            errorFlag = 9;
        g = common[0];

        // Additional error checking including preventing two nodes connecting
        // multiple times
        if (c == d || e == f)
            errorFlag = 6; // can simply be triangle edge sharing pair (not an error)
        if (vContains(networkB.nodes[w].netCnxs, x))
            errorFlag = 7;
        if (vContains(networkB.nodes[y].netCnxs, x))
            errorFlag = 7;
        if (vContains(networkB.nodes[y].auxCnxs, x))
            errorFlag = 7;

        if (errorFlag != 0) {
            return false;
        }

        // check move will not violate dual connectivity limits
        if (networkB.nodes[u].netCnxs.n == minBCnxs ||
            networkB.nodes[v].netCnxs.n == minBCnxs ||
            networkB.nodes[w].netCnxs.n == maxBCnxs ||
            networkB.nodes[x].netCnxs.n == maxBCnxs)
            return false;
        else {
            switchIdsA = VecF<int>(7);
            switchIdsB = VecF<int>(5);
            switchIdsA[0] = a;
            switchIdsA[1] = b;
            switchIdsA[2] = c;
            switchIdsA[3] = d;
            switchIdsA[4] = e;
            switchIdsA[5] = f;
            switchIdsA[6] = g;
            switchIdsB[0] = u;
            switchIdsB[1] = v;
            switchIdsB[2] = w;
            switchIdsB[3] = x;
            switchIdsB[4] = y;
            return true;
        }
    } else {
        throw std::runtime_error("Not yet implemented");
    }
}

const int MIN_COORDINATION_NUMBER = 2;
const int NUM_MIX_IDS_A = 6;
const int NUM_MIX_IDS_B = 7;

/**
 * @brief Generate all ids of nodes in lattices A and B needed for mix move, for
 * all coordinations >=2
 * @param cnxType Connection type of nodes to be mixed
 * @param mixIdsA Vector of node ids in lattice A
 * @param mixIdsB Vector of node ids in lattice B
 * @param a Node id in lattice A
 * @param b Node id in lattice A
 * @param u Node id in lattice B
 * @param v Node id in lattice B
 * @return True if move is possible, false otherwise
 */
bool LinkedNetwork::generateMixIds(int cnxType, VecF<int> &mixIdsA,
                                   VecF<int> &mixIdsB, int a, int b, int u,
                                   int v) {
    if (cnxType <
        MIN_COORDINATION_NUMBER) { // cannot mix if first node is 2 coordinate
        return false;
    }
    /* Get required node ids in lattice and dual
     * a,b,c,d,e,f are nodes in lattice A
     * u,v,w,x,y,z,xx are nodes in lattice B
     *  E      F            V
     *   \    /          /  |  \
     * --A---B--    XX--X   |   Z
     *  /     \         W   |   Y
     * C       D         \  |  /
     *                      U
     */

    int errorFlag = 0;
    int c, d, e, f;
    int w, x, y, z, xx;

    VecR<int> common, common1;
    // c is defined as sharing a,u but not being b
    c = findCommonConnection(a, u, b);
    // d is defined as sharing b,u but not being a
    d = findCommonConnection(b, u, a);
    // e is defined as sharing a,v but not being b
    e = findCommonConnection(a, v, b);
    // f is defined as sharing b,v but not being a
    f = findCommonConnection(b, v, a);

    // w is defines as sharing a,c but not being u
    w = findAssociatedNodeAA(a, c, u);
    // x is defines as sharing a,e but not being v
    x = findAssociatedNodeAA(a, e, v);
    // y is defines as sharing b,d but not being u
    y = findAssociatedNodeAA(b, d, u);
    // z is defines as sharing b,f but not being v
    z = findAssociatedNodeAA(b, f, v);

    // Find next node connected to x by inspecting around a
    int nCnxs = networkA.nodes[a].dualCnxs.n;
    for (int i = 0; i < nCnxs; ++i) {
        if (networkA.nodes[a].dualCnxs[i] == u) {
            int j = (i + 1) % nCnxs;
            int k = (i - 1 + nCnxs) % nCnxs;
            if (networkA.nodes[a].dualCnxs[j] == v)
                xx = networkA.nodes[a].dualCnxs[(j + 2) % nCnxs];
            else if (networkA.nodes[a].dualCnxs[k] == v)
                xx = networkA.nodes[a].dualCnxs[(k - 2 + nCnxs) % nCnxs];
            else
                std::cout << "Error in xx detection" << std::endl;
        }
    }

    // Additional error checking
    if (c == d || e == f)
        errorFlag = 6; // can simply be triangle edge sharing pair (not an error)
    // Prevent rings having only two or fewer neighbours
    VecR<int> vCnxs = vUnique(networkB.nodes[v].netCnxs);
    if (vCnxs.n <= 3)
        errorFlag = 10;
    VecR<int> uCnxs = vUnique(networkB.nodes[u].netCnxs);
    uCnxs.delValue(v);
    if (uCnxs.n <= 2) {
        for (int i = 0; i < uCnxs.n; ++i) {
            if (uCnxs[i] == x) {
                errorFlag = 10;
                break;
            }
        }
    }
    if (errorFlag != 0)
        return false;
    // check move will not violate connectivity limits
    if (networkB.nodes[v].netCnxs.n == minBCnxs ||
        networkB.nodes[x].netCnxs.n == maxBCnxs ||
        networkA.nodes[a].netCnxs.n == minACnxs ||
        networkA.nodes[b].netCnxs.n == maxACnxs) {
        return false;
    }
    mixIdsA = VecF<int>(NUM_MIX_IDS_A);
    mixIdsA[0] = a;
    mixIdsA[1] = b;
    mixIdsA[2] = c;
    mixIdsA[3] = d;
    mixIdsA[4] = e;
    mixIdsA[5] = f;

    mixIdsB = VecF<int>(NUM_MIX_IDS_B);
    mixIdsB[0] = u;
    mixIdsB[1] = v;
    mixIdsB[2] = w;
    mixIdsB[3] = x;
    mixIdsB[4] = y;
    mixIdsB[5] = z;
    mixIdsB[6] = xx;
    return true;
}

int LinkedNetwork::findCommonConnection(int node1, int node2, int excludeNode) {

    // Find node that shares node1 and node2 but is not excludeNode
    int commonConnection = -1;
    VecR<int> commonConnections = vCommonValues(networkA.nodes[node1].netCnxs, networkB.nodes[node2].dualCnxs);
    commonConnections.delValue(excludeNode);
    if (commonConnections.n == 1)
        commonConnection = commonConnections[0];
    else { // rare high temperature occurrence as a result of 2-cnd nodes giving
           // ring inside ring
        VecR<int> common1(0, commonConnections.n);
        int nCnxs = networkA.nodes[node1].netCnxs.n;
        for (int i = 0; i < nCnxs; ++i) {
            if (networkA.nodes[node1].netCnxs[i] == excludeNode) {
                int l = networkA.nodes[node1].netCnxs[(i + 1) % nCnxs];
                int r = networkA.nodes[node1].netCnxs[(i - 1 + nCnxs) % nCnxs];
                for (int j = 0; j < commonConnections.n; ++j) {
                    if (commonConnections[j] == l)
                        common1.addValue(commonConnections[j]);
                    else if (commonConnections[j] == r)
                        common1.addValue(commonConnections[j]);
                }
                break;
            }
        }
        if (common1.n == 1)
            commonConnection = common1[0];
        else { // even rarer case
            int nCnxs = networkB.nodes[node2].dualCnxs.n;
            for (int i = 0; i < nCnxs; ++i) {
                if (networkB.nodes[node2].dualCnxs[i] == excludeNode) {
                    int j;
                    if (networkB.nodes[node2].dualCnxs[(i + 1) % nCnxs] == node1)
                        j = (i + 2) % nCnxs;
                    else if (networkB.nodes[node2].dualCnxs[(i + nCnxs - 1) % nCnxs] == node1)
                        j = (i + nCnxs - 2) % nCnxs;
                    if (vContains(common1, networkB.nodes[node2].dualCnxs[j])) {
                        commonConnection = networkB.nodes[node2].dualCnxs[j];
                        break;
                    }
                }
            }
        }
    }

    if (commonConnection == -1) {
        wrapCoordinates();
        syncCoordinates();
        write("debug");
        std::ostringstream oss;
        oss << "Could not find associated node for nodes " << node1 << " and " << node2
            << " with exclude node " << excludeNode;
        throw std::runtime_error(oss.str());
    }

    return commonConnection;
}

int LinkedNetwork::findAssociatedNodeAA(int idA, int idB, int idDel) {

    // Find node that shares idA and idB but is not idDel
    int associated = -1;
    VecR<int> common;
    common =
        vCommonValues(networkA.nodes[idA].dualCnxs, networkA.nodes[idB].dualCnxs);
    common.delValue(idDel);
    if (common.n == 1)
        associated = common[0];
    else { // rare case with periodic interactions
        VecR<int> common1(0, common.n);
        for (int i = 0; i < common.n; ++i) {
            if (vContains(networkB.nodes[common[i]].netCnxs, idDel))
                common1.addValue(common[i]);
        }
        if (common1.n == 1)
            associated = common1[0];
        else { // rare case of large ring surrounding group
            // check a,b adjacent on ring
            VecR<int> common2(0, common1.n);
            for (int i = 0; i < common1.n; ++i) {
                VecR<int> ring = networkB.nodes[common[i]].dualCnxs;
                for (int j = 0; j < ring.n; ++j) {
                    int k = (j + 1) % ring.n;
                    if (ring[j] == idA && ring[k] == idB) {
                        common2.addValue(common1[i]);
                        break;
                    } else if (ring[j] == idB && ring[k] == idA) {
                        common2.addValue(common1[i]);
                        break;
                    }
                }
            }
            if (common2.n == 1)
                associated = common2[0];
        }
    }

    if (associated == -1) {
        wrapCoordinates();
        syncCoordinates();
        write("debug");
        std::cout << idA << " " << idB << " " << idDel << " " << common
                  << std::endl;
        std::cout << "ERROR IN ASSOCIATED NODE" << std::endl;
    }

    return associated;
}

// Switch connectivities in lattice between 2x3 coordinate nodes
void LinkedNetwork::switchCnx33(VecF<int> switchIdsA, VecF<int> switchIdsB,
                                VecF<int> switchIdsT) {
    // unpck parameters
    int a = switchIdsA[0];
    int b = switchIdsA[1];
    int c = switchIdsA[2];
    int d = switchIdsA[3];
    int e = switchIdsA[4];
    int f = switchIdsA[5];
    int u = switchIdsB[0];
    int v = switchIdsB[1];
    int w = switchIdsB[2];
    int x = switchIdsB[3];

    int alpha = switchIdsT[6];
    int beta = switchIdsT[7];
    int gamma = switchIdsT[8];
    int delta = switchIdsT[9];
    int eta = switchIdsT[10];
    // Apply changes to descriptors due to breaking connections
    // For network A node distribution and edge distribution will remain unchanged

    // For network B node and edge distribution will change
    int nu = networkB.nodes[u].netCnxs.n;
    int nv = networkB.nodes[v].netCnxs.n;
    int nw = networkB.nodes[w].netCnxs.n;
    int nx = networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[u].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if (id != v && id != w && id != x)
            --networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if (id != u && id != w && id != x)
            --networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if (id != u && id != v && id != x)
            --networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if (id != u && id != v && id != w)
            --networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }

    // A-A connectivities
    // swap a->e/d, b->d/e, d->b/a, e->a/b
    networkA.nodes[a].netCnxs.swapValue(e, d);
    networkA.nodes[b].netCnxs.swapValue(d, e);
    networkA.nodes[d].netCnxs.swapValue(b, a);
    networkA.nodes[e].netCnxs.swapValue(a, b);

    // A-B connectvities
    // swap a->v/x, b->u/w
    networkA.nodes[a].dualCnxs.swapValue(v, x);
    networkA.nodes[b].dualCnxs.swapValue(u, w);

    // B-B connectivities
    // have to account for the fact that two nodes may connect multiple times
    // break insert w:u-(x)-v, x:u-(w)-v

    networkB.nodes[u].netCnxs.swapValue(v, -1, w, x);
    networkB.nodes[v].netCnxs.swapValue(u, -1, w, x);
    networkB.nodes[u].netCnxs.delValue(-1);
    networkB.nodes[v].netCnxs.delValue(-1);
    if (vAdjCount(networkB.nodes[w].netCnxs, u, v) > 1) {
        VecR<int> eCnxs = networkA.nodes[e].dualCnxs;
        eCnxs.delValue(v);
        eCnxs.delValue(w);
        int ww = eCnxs[0];
        networkB.nodes[w].netCnxs.swapValue(v, -1, ww, u);
        networkB.nodes[w].netCnxs.insertValue(x, -1, u);
        networkB.nodes[w].netCnxs.swapValue(-1, v);
    } else {
        networkB.nodes[w].netCnxs.insertValue(x, u, v);
    }
    if (vAdjCount(networkB.nodes[x].netCnxs, u, v) > 1) {
        VecR<int> fCnxs = networkA.nodes[f].dualCnxs;
        fCnxs.delValue(v);
        fCnxs.delValue(x);
        int xx = fCnxs[0];
        networkB.nodes[x].netCnxs.swapValue(v, -1, xx, u);
        networkB.nodes[x].netCnxs.insertValue(w, -1, u);
        networkB.nodes[x].netCnxs.swapValue(-1, v);
    } else {
        networkB.nodes[x].netCnxs.insertValue(w, u, v);
    }

    // B-A connectivities
    // break u->b, v->a, insert w:a-(b)-e, z:b-(a)-d
    networkB.nodes[u].dualCnxs.delValue(b);
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b, a, e);
    networkB.nodes[x].dualCnxs.insertValue(a, b, d);

    // Apply changes to descriptors due to making connections
    // Network B
    nu = networkB.nodes[u].netCnxs.n;
    nv = networkB.nodes[v].netCnxs.n;
    nw = networkB.nodes[w].netCnxs.n;
    nx = networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[u].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if (id != v && id != w && id != x)
            ++networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if (id != u && id != w && id != x)
            ++networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if (id != u && id != v && id != x)
            ++networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if (id != u && id != v && id != w)
            ++networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }

    networkT.nodes[a].netCnxs.swapValue(gamma, delta);
    networkT.nodes[b].netCnxs.swapValue(delta, gamma);
    networkT.nodes[delta].netCnxs.swapValue(b, a);
    networkT.nodes[gamma].netCnxs.swapValue(a, b);
    networkT.nodes[gamma].netCnxs.swapValue(beta, eta);
    networkT.nodes[beta].netCnxs.swapValue(gamma, delta);
    networkT.nodes[delta].netCnxs.swapValue(eta, beta);
    networkT.nodes[eta].netCnxs.swapValue(delta, gamma);
    networkT.nodes[a].dualCnxs.swapValue(v, x);
    networkT.nodes[b].dualCnxs.swapValue(u, w);
}

// Switch connectivities in lattice between 2x4 coordinate nodes
void LinkedNetwork::switchCnx44(VecF<int> switchIdsA, VecF<int> switchIdsB) {

    // unpck parameters
    int a = switchIdsA[0];
    int b = switchIdsA[1];
    int c = switchIdsA[2];
    int d = switchIdsA[3];
    int e = switchIdsA[4];
    int f = switchIdsA[5];
    int g = switchIdsA[6];
    int h = switchIdsA[7];
    int u = switchIdsB[0];
    int v = switchIdsB[1];
    int w = switchIdsB[2];
    int x = switchIdsB[3];
    int y = switchIdsB[4];
    int z = switchIdsB[5];

    // Apply changes to descriptors due to breaking connections
    // For network A node distribution and edge distribution will remain unchanged

    // For network B node and edge distribution will change
    int nu = networkB.nodes[u].netCnxs.n;
    int nv = networkB.nodes[v].netCnxs.n;
    int nw = networkB.nodes[w].netCnxs.n;
    int nx = networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[u].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if (id != v && id != w && id != x)
            --networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if (id != u && id != w && id != x)
            --networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if (id != u && id != v && id != x)
            --networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if (id != u && id != v && id != w)
            --networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }

    // A-A connectivities
    // break a->e, b->d, insert a:b-(d)-c b:a-(e)-f, swap d->b/a, e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.delValue(d);
    networkA.nodes[a].netCnxs.insertValue(d, b, c);
    networkA.nodes[b].netCnxs.insertValue(e, a, f);
    networkA.nodes[d].netCnxs.swapValue(b, a);
    networkA.nodes[e].netCnxs.swapValue(a, b);

    // A-B connectvities
    // swap a->v/x, b->u/w
    networkA.nodes[a].dualCnxs.swapValue(v, x);
    networkA.nodes[b].dualCnxs.swapValue(u, w);

    // B-B connectivities
    // have to account for the fact that two nodes may connect multiple times
    // break u:x-(v)-y, v:w-(u)-z, insert w:y-(x)-v, x:u-(w)-z
    networkB.nodes[u].netCnxs.swapValue(v, -1, x, y);
    networkB.nodes[v].netCnxs.swapValue(u, -1, w, z);
    networkB.nodes[u].netCnxs.delValue(-1);
    networkB.nodes[v].netCnxs.delValue(-1);
    networkB.nodes[w].netCnxs.insertValue(x, y, v);
    networkB.nodes[x].netCnxs.insertValue(w, u, z);

    // B-A connectivities
    // break u->b, v->a, insert w:a-(b)-e, x:b-(a)-d
    networkB.nodes[u].dualCnxs.delValue(b);
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b, a, e);
    networkB.nodes[x].dualCnxs.insertValue(a, b, d);

    // B-B secondary connectivities
    // break v-y, u-z, swap y->v/x, z->u/w, add x-y, w-z
    networkB.nodes[v].auxCnxs.delValue(y);
    networkB.nodes[u].auxCnxs.delValue(z);
    networkB.nodes[y].auxCnxs.swapValue(v, x);
    networkB.nodes[z].auxCnxs.swapValue(u, w);
    networkB.nodes[x].auxCnxs.addValue(y);
    networkB.nodes[w].auxCnxs.addValue(z);

    // Apply changes to descriptors due to making connections
    // Network B
    nu = networkB.nodes[u].netCnxs.n;
    nv = networkB.nodes[v].netCnxs.n;
    nw = networkB.nodes[w].netCnxs.n;
    nx = networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[u].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if (id != v && id != w && id != x)
            ++networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if (id != u && id != w && id != x)
            ++networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if (id != u && id != v && id != x)
            ++networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if (id != u && id != v && id != w)
            ++networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }
}

// Switch connectivities in lattice between 4 and 3 coordinate nodes
void LinkedNetwork::switchCnx43(VecF<int> switchIdsA, VecF<int> switchIdsB) {

    // unpck parameters
    int a = switchIdsA[0];
    int b = switchIdsA[1];
    int c = switchIdsA[2];
    int d = switchIdsA[3];
    int e = switchIdsA[4];
    int f = switchIdsA[5];
    int g = switchIdsA[6];
    int u = switchIdsB[0];
    int v = switchIdsB[1];
    int w = switchIdsB[2];
    int x = switchIdsB[3];
    int y = switchIdsB[4];

    // Apply changes to descriptors due to breaking connections
    // For network A node distribution will remain unchanged but edge distribution
    // will change
    int na = networkA.nodes[a].netCnxs.n;
    int nb = networkA.nodes[b].netCnxs.n;
    int nd = networkA.nodes[d].netCnxs.n;
    int ne = networkA.nodes[e].netCnxs.n;
    --networkA.edgeDistribution[na][ne];
    --networkA.edgeDistribution[nb][nd];
    --networkA.edgeDistribution[ne][na];
    --networkA.edgeDistribution[nd][nb];

    // For network B node and edge distribution will change
    int nu = networkB.nodes[u].netCnxs.n;
    int nv = networkB.nodes[v].netCnxs.n;
    int nw = networkB.nodes[w].netCnxs.n;
    int nx = networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[u].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if (id != v && id != w && id != x)
            --networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if (id != u && id != w && id != x)
            --networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if (id != u && id != v && id != x)
            --networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if (id != u && id != v && id != w)
            --networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }

    // A-A connectivities
    // break a->e, insert a:b-(d)-c, swap b->d/e, d->b/a, e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.swapValue(d, e);
    networkA.nodes[a].netCnxs.insertValue(d, b, c);
    networkA.nodes[d].netCnxs.swapValue(b, a);
    networkA.nodes[e].netCnxs.swapValue(a, b);

    // A-B connectvities
    // swap a->v/x, b->u/w
    networkA.nodes[a].dualCnxs.swapValue(v, x);
    networkA.nodes[b].dualCnxs.swapValue(u, w);

    // B-B connectivities
    // have to account for the fact that two nodes may connect multiple times
    // break u:x-(v)-y, v:w-(u)-x, insert w:y-(x)-v, x:u-(w)-v
    networkB.nodes[u].netCnxs.swapValue(v, -1, x, y);
    networkB.nodes[v].netCnxs.swapValue(u, -1, w, x);
    networkB.nodes[u].netCnxs.delValue(-1);
    networkB.nodes[v].netCnxs.delValue(-1);
    networkB.nodes[w].netCnxs.insertValue(x, y, v);
    networkB.nodes[x].netCnxs.insertValue(w, u, v);

    // B-A connectivities
    // break u->b, v->a, insert w:a-(b)-e, x:b-(a)-d
    networkB.nodes[u].dualCnxs.delValue(b);
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b, a, e);
    networkB.nodes[x].dualCnxs.insertValue(a, b, d);

    // B-B secondary connectivities
    // break v-y, swap y->v/x, add x-y
    networkB.nodes[v].auxCnxs.delValue(y);
    networkB.nodes[y].auxCnxs.swapValue(v, x);
    networkB.nodes[x].auxCnxs.addValue(y);

    // Apply changes to descriptors due to making connections
    // Network A
    ++networkA.edgeDistribution[na][nd];
    ++networkA.edgeDistribution[nb][ne];
    ++networkA.edgeDistribution[nd][na];
    ++networkA.edgeDistribution[ne][nb];
    // Network B
    nu = networkB.nodes[u].netCnxs.n;
    nv = networkB.nodes[v].netCnxs.n;
    nw = networkB.nodes[w].netCnxs.n;
    nx = networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[u].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if (id != v && id != w && id != x)
            ++networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if (id != u && id != w && id != x)
            ++networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if (id != u && id != v && id != x)
            ++networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if (id != u && id != v && id != w)
            ++networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }
}

// Mix connectivities to exchange 3/4 coordination nodes
void LinkedNetwork::mixCnx34(VecF<int> mixIdsA, VecF<int> mixIdsB) {
    // unpck parameters
    int a = mixIdsA[0];
    int b = mixIdsA[1];
    int c = mixIdsA[2];
    int d = mixIdsA[3];
    int e = mixIdsA[4];
    int f = mixIdsA[5];
    int g = mixIdsA[6];
    int u = mixIdsB[0];
    int v = mixIdsB[1];
    int w = mixIdsB[2];
    int x = mixIdsB[3];
    int y = mixIdsB[4];

    // Apply changes to descriptors due to breaking connections
    // For network A node distribution will remain unchanged but edge distribution
    // will change
    int na = networkA.nodes[a].netCnxs.n;
    int nb = networkA.nodes[b].netCnxs.n;
    for (int i = 0; i < na; ++i) {
        int id = networkA.nodes[a].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[na][nCnx];
        if (id != b)
            --networkA.edgeDistribution[nCnx][na]; // prevent double counting
    }
    for (int i = 0; i < nb; ++i) {
        int id = networkA.nodes[b].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[nb][nCnx];
        if (id != a)
            --networkA.edgeDistribution[nCnx][nb]; // prevent double counting
    }

    // For network B node and edge distribution will change
    int nv = networkB.nodes[v].netCnxs.n;
    int nw = networkB.nodes[w].netCnxs.n;
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if (id != w)
            --networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if (id != v)
            --networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }

    // A-A connectivities
    // break a->e, insert b:a-(e)-f, swap e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.insertValue(e, a, f);
    networkA.nodes[e].netCnxs.swapValue(a, b);

    // A-B connectvities
    // break a->v, insert b:u-(w)-v
    networkA.nodes[a].dualCnxs.delValue(v);
    networkA.nodes[b].dualCnxs.insertValue(w, u, v);

    // B-B connectivities
    // break v->u, insert w:y-(u)-v, swap u->v/w
    networkB.nodes[v].netCnxs.delValue(u);
    networkB.nodes[w].netCnxs.insertValue(u, y, v);
    networkB.nodes[u].netCnxs.swapValue(v, w);

    // B-A connectivities
    // break v->a, insert w:a-(b)-e
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[w].dualCnxs.insertValue(b, a, e);

    // B-B secondary connectivities
    // break y-v, u-z, swap u->w/v, v->y/u, w->u/x, add x-w
    networkB.nodes[y].auxCnxs.delValue(v);
    networkB.nodes[u].auxCnxs.swapValue(w, v);
    networkB.nodes[w].auxCnxs.swapValue(u, x);
    networkB.nodes[v].auxCnxs.swapValue(y, u);
    networkB.nodes[x].auxCnxs.addValue(w);

    // Apply changes to descriptors due to making connections
    // Network A
    na = networkA.nodes[a].netCnxs.n;
    nb = networkA.nodes[b].netCnxs.n;
    for (int i = 0; i < na; ++i) {
        int id = networkA.nodes[a].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[na][nCnx];
        if (id != b)
            ++networkA.edgeDistribution[nCnx][na]; // prevent double counting
    }
    for (int i = 0; i < nb; ++i) {
        int id = networkA.nodes[b].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[nb][nCnx];
        if (id != a)
            ++networkA.edgeDistribution[nCnx][nb]; // prevent double counting
    }
    // Network B
    nv = networkB.nodes[v].netCnxs.n;
    nw = networkB.nodes[w].netCnxs.n;
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if (id != w)
            ++networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[w].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if (id != v)
            ++networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
}

// Mix connectivities of adjacent nodes to give -1/+1 in coordinations
void LinkedNetwork::mixCnx(VecF<int> mixIdsA, VecF<int> mixIdsB) {

    // unpack parameters
    /* a,b,c,d,e,f are nodes in lattice A
     * u,v,w,x,y,z are nodes in lattice B
     *  E      F            V
     *   \    /          /  |  \
     * --A---B--    XX--X   |   Z
     *  /     \         W   |   Y
     * C       D         \  |  /
     *                      U
     *
     *       E F               V
     *       |/              /  \
     * --A---B--        XX--X   Z
     *  /     \         W   |   Y
     * C       D         \  |  /
     *                      U
     */
    int a, b, c, d, e, f;
    int u, v, w, x, y, z, xx;
    a = mixIdsA[0];
    b = mixIdsA[1];
    c = mixIdsA[2];
    d = mixIdsA[3];
    e = mixIdsA[4];
    f = mixIdsA[5];
    u = mixIdsB[0];
    v = mixIdsB[1];
    w = mixIdsB[2];
    x = mixIdsB[3];
    y = mixIdsB[4];
    z = mixIdsB[5];
    xx = mixIdsB[6];

    // Apply changes to descriptors due to breaking connections
    // For network A node distribution and edge distribution will change as a
    // result of a,b mixing
    int na, nb;
    na = networkA.nodes[a].netCnxs.n;
    nb = networkA.nodes[b].netCnxs.n;
    --networkA.nodeDistribution[na];
    --networkA.nodeDistribution[nb];
    for (int i = 0; i < na; ++i) {
        int id = networkA.nodes[a].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[na][nCnx];
        if (id != b)
            --networkA.edgeDistribution[nCnx][na]; // prevent double counting
    }
    for (int i = 0; i < nb; ++i) {
        int id = networkA.nodes[b].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        --networkA.edgeDistribution[nb][nCnx];
        if (id != a)
            --networkA.edgeDistribution[nCnx][nb]; // prevent double counting
    }

    // For network B node and edge distribution will change as a result of v,x
    // mixing
    int nv, nx;
    nv = networkB.nodes[v].netCnxs.n;
    nx = networkB.nodes[x].netCnxs.n;
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nx];
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if (id != x)
            --networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if (id != v)
            --networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }

    // A-A connectivities
    // break a->e, insert b:a-(e)-f, swap e->a/b
    networkA.nodes[a].netCnxs.delValue(e);
    networkA.nodes[b].netCnxs.insertValue(e, a, f);
    networkA.nodes[e].netCnxs.swapValue(a, b);

    // A-B connectvities
    // break a->v, insert b:u-(x)-v
    networkA.nodes[a].dualCnxs.delValue(v);
    networkA.nodes[b].dualCnxs.insertValue(x, u, v);

    // B-B connectivities
    // break v->u, insert x:xx-(u)-v, swap u->v/x
    networkB.nodes[v].netCnxs.swapValue(
        u, -1, x, z); // swap and delete as can occur multiple times if 2-cnd node
    networkB.nodes[v].netCnxs.delValue(-1);
    networkB.nodes[x].netCnxs.insertValue(u, xx, v);
    networkB.nodes[u].netCnxs.swapValue(v, x, w, y);

    // B-A connectivities
    // break v->a, insert x:a-(b)-e
    networkB.nodes[v].dualCnxs.delValue(a);
    networkB.nodes[x].dualCnxs.insertValue(b, a, e);

    // B-B secondary connectivities
    // DONT NEED FOR THIS GENERAL MIX MOVE

    // Apply changes to descriptors due to making connections
    // Network A
    na = networkA.nodes[a].netCnxs.n;
    nb = networkA.nodes[b].netCnxs.n;
    ++networkA.nodeDistribution[na];
    ++networkA.nodeDistribution[nb];
    for (int i = 0; i < na; ++i) {
        int id = networkA.nodes[a].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[na][nCnx];
        if (id != b)
            ++networkA.edgeDistribution[nCnx][na]; // prevent double counting
    }
    for (int i = 0; i < nb; ++i) {
        int id = networkA.nodes[b].netCnxs[i];
        int nCnx = networkA.nodes[id].netCnxs.n;
        ++networkA.edgeDistribution[nb][nCnx];
        if (id != a)
            ++networkA.edgeDistribution[nCnx][nb]; // prevent double counting
    }
    // Network B
    nv = networkB.nodes[v].netCnxs.n;
    nx = networkB.nodes[x].netCnxs.n;
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nx];
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[v].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if (id != x)
            ++networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[x].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if (id != v)
            ++networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }
}

// Check mix move did not introduce any edges which form part of three rings
bool LinkedNetwork::checkThreeRingEdges(int id) {

    bool edgeCheck = true;
    int nCnxs = networkB.nodes[id].dualCnxs.n;
    for (int i = 0; i < nCnxs; ++i) {
        int j = networkB.nodes[id].dualCnxs[i];
        int k = networkB.nodes[id].dualCnxs[(i + 1) % nCnxs];
        VecR<int> common =
            vCommonValues(networkA.nodes[j].dualCnxs, networkA.nodes[k].dualCnxs);
        if (common.n > 2) {
            edgeCheck = false;
            break;
        }
    }
    return edgeCheck;
}

// Rearrange nodes after connection switch to maintain convexity
bool LinkedNetwork::convexRearrangement(int cnxType, VecF<int> switchIdsA,
                                        VecF<int> switchIdsB) {

    bool convex;

    if (cnxType == 33 || cnxType == 34 || cnxType == 44) {
        // Unpack nodes
        int a, b, c, d, e, f;
        a = switchIdsA[0];
        b = switchIdsA[1];
        c = switchIdsA[2];
        d = switchIdsA[3];
        e = switchIdsA[4];
        f = switchIdsA[5];

        // Maintains which nodes are convex
        VecF<bool> convexNodes;
        if (cnxType == 33)
            convexNodes = VecF<bool>(6);
        else if (cnxType == 34)
            convexNodes = VecF<bool>(7);
        else if (cnxType == 44)
            convexNodes = VecF<bool>(8);

        // Initial guess places a at the centreCoords of cd, b at the centreCoords of ef
        VecF<double> va(2), vb(2), vc(2), vd(2), ve(2), vf(2);
        va[0] = crds[2 * a];
        va[1] = crds[2 * a + 1];
        vb[0] = crds[2 * b];
        vb[1] = crds[2 * b + 1];
        vc[0] = crds[2 * c];
        vc[1] = crds[2 * c + 1];
        vd[0] = crds[2 * d];
        vd[1] = crds[2 * d + 1];
        ve[0] = crds[2 * e];
        ve[1] = crds[2 * e + 1];
        vf[0] = crds[2 * f];
        vf[1] = crds[2 * f + 1];
        VecF<double> vce(2), vdf(2), vcd(2), vef(2);
        vce = ve - vc;
        vdf = vf - vd;
        vcd = vd - vc;
        vef = vf - ve;
        vce[0] -= networkA.dimensions[0] * nearbyint(vce[0] * networkA.rpb[0]);
        vce[1] -= networkA.dimensions[1] * nearbyint(vce[1] * networkA.rpb[1]);
        vdf[0] -= networkA.dimensions[0] * nearbyint(vdf[0] * networkA.rpb[0]);
        vdf[1] -= networkA.dimensions[1] * nearbyint(vdf[1] * networkA.rpb[1]);
        vcd[0] -= networkA.dimensions[0] * nearbyint(vcd[0] * networkA.rpb[0]);
        vcd[1] -= networkA.dimensions[1] * nearbyint(vcd[1] * networkA.rpb[1]);
        vef[0] -= networkA.dimensions[0] * nearbyint(vef[0] * networkA.rpb[0]);
        vef[1] -= networkA.dimensions[1] * nearbyint(vef[1] * networkA.rpb[1]);
        va = vd - vcd / 2.0 + vdf / 10.0;
        vb = ve + vef / 2.0 - vce / 10.0;
        va[0] -= networkA.dimensions[0] * nearbyint(va[0] * networkA.rpb[0]);
        va[1] -= networkA.dimensions[1] * nearbyint(va[1] * networkA.rpb[1]);
        vb[0] -= networkA.dimensions[0] * nearbyint(vb[0] * networkA.rpb[0]);
        vb[1] -= networkA.dimensions[1] * nearbyint(vb[1] * networkA.rpb[1]);
        crds[2 * a] = va[0];
        crds[2 * a + 1] = va[1];
        crds[2 * b] = vb[0];
        crds[2 * b + 1] = vb[1];
        for (int i = 0; i < switchIdsA.n; ++i)
            convexNodes[i] = checkConvexity(switchIdsA[i]);
        convex = (convexNodes == true);

        // Guess move a,b towards each other
        VecF<double> vab(2);
        if (!convex) {
            vab = (vb - va) * 0.5 * 0.1;
            vab[0] -= networkA.dimensions[0] * nearbyint(vab[0] * networkA.rpb[0]);
            vab[1] -= networkA.dimensions[1] * nearbyint(vab[1] * networkA.rpb[1]);
            for (int i = 0; i < 9; ++i) {
                va += vab;
                vb -= vab;
                va[0] -= networkA.dimensions[0] * nearbyint(va[0] * networkA.rpb[0]);
                va[1] -= networkA.dimensions[1] * nearbyint(va[1] * networkA.rpb[1]);
                vb[0] -= networkA.dimensions[0] * nearbyint(vb[0] * networkA.rpb[0]);
                vb[1] -= networkA.dimensions[1] * nearbyint(vb[1] * networkA.rpb[1]);
                crds[2 * a] = va[0];
                crds[2 * a + 1] = va[1];
                crds[2 * b] = vb[0];
                crds[2 * b + 1] = vb[1];
                for (int i = 0; i < switchIdsA.n; ++i)
                    convexNodes[i] = checkConvexity(switchIdsA[i]);
                convex = (convexNodes == true);
                if (convex)
                    break;
            }
        }
    } else
        throw(std::string("Not yet implemented!"));

    return convex;
}

// Single monte carlo switching move
VecF<int> LinkedNetwork::monteCarloSwitchMoveLAMMPS(double &SimpleGrapheneEnergy, double &TersoffGrapheneEnergy,
                                                    double &TriangleRaftEnergy, double &BilayerEnergy, double &BNEnergy,
                                                    LoggerPtr logger) {

    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */

    // Select valid random connection - that will not violate connection limits
    int a;
    int b;
    int u;
    int v;
    VecF<int> switchIdsA;
    VecF<int> switchIdsB;
    VecF<int> switchIdsT;
    bool foundValidMove = false;
    int cnxType;
    logger->debug("Finding move...");
    for (int i = 0; i < networkA.nodes.n * networkA.nodes.n; ++i) {
        if (mcWeighting == "weighted")
            cnxType = pickRandomConnection(a, b, u, v, mtGen, SelectionType::EXPONENTIAL_DECAY, logger);
        else
            cnxType = pickRandomConnection(a, b, u, v, mtGen, SelectionType::RANDOM, logger);
        foundValidMove = generateSwitchIDs(cnxType, switchIdsA, switchIdsB, switchIdsT, a, b, u, v, logger);
        if (foundValidMove)
            break;
    }
    if (!foundValidMove) {
        logger->error("Cannot find any valid switch moves");
        throw std::runtime_error("Cannot find any valid switch moves");
    }

    // Save current state
    logger->debug("Saving energies...");
    double saveEnergySimpleGraphene = SimpleGrapheneEnergy;
    double saveEnergyTersoffGraphene = TersoffGrapheneEnergy;
    double saveEnergyTriangleRaft = TriangleRaftEnergy;
    double saveEnergyBilayer = BilayerEnergy;
    double saveEnergyBN = BNEnergy;
    if (isSimpleGrapheneEnabled) {
        if (abs(SimpleGrapheneEnergy - SimpleGraphene.globalPotentialEnergy()) > 0.001) {
            logger->debug("Saved simpleGraphene energy: {} vs calculated: {}",
                          SimpleGrapheneEnergy,
                          SimpleGraphene.globalPotentialEnergy());
        }
    }
    if (isTriangleRaftEnabled) {
        if (abs(TriangleRaftEnergy - Triangle_Raft.globalPotentialEnergy()) >
            0.001) {
            logger->debug("Saved triangleRaft energy: {} vs calculated: {}",
                          TriangleRaftEnergy, Triangle_Raft.globalPotentialEnergy());
        } else {
            logger->debug("Triangle Raft : {}",
                          Triangle_Raft.globalPotentialEnergy());
        }
    }
    logger->debug("Saving coordinates...");
    VecF<double> saveCrds = crds;
    double *saveCrdsSimpleGraphene;
    double *saveCrdsTersoffGraphene;
    double *saveCrdsTriangleRaft;
    double *saveCrdsBilayer;
    double *saveCrdsBN;

    if (isSimpleGrapheneEnabled) {
        logger->debug("Saving SimpleGraphene coordinates");
        saveCrdsSimpleGraphene = SimpleGraphene.fetchCrds(2);
    }
    if (isTersoffGrapheneEnabled) {
        logger->debug("Saving TersoffGraphene coordinates");
        saveCrdsTersoffGraphene = TersoffGraphene.fetchCrds(2);
    }
    if (isTriangleRaftEnabled) {
        logger->debug("Saving TriangleRaft coordinates");
        saveCrdsTriangleRaft = Triangle_Raft.fetchCrds(2);
    }
    if (isBilayerEnabled) {
        logger->debug("Saving Bilayer coordinates");
        saveCrdsBilayer = Bilayer.fetchCrds(3);
    }
    if (isBNEnabled) {
        logger->debug("Saving BN coordinates");
        saveCrdsBN = BN.fetchCrds(2);
    }
    VecF<int> saveNodeDistA = networkA.nodeDistribution;
    VecF<int> saveNodeDistB = networkB.nodeDistribution;

    VecF<VecF<int>> saveEdgeDistA = networkA.edgeDistribution;
    VecF<VecF<int>> saveEdgeDistB = networkB.edgeDistribution;

    VecF<Node> saveNodesA(switchIdsA.n);
    VecF<Node> saveNodesB(switchIdsB.n);
    VecF<Node> saveNodesT(switchIdsT.n);
    for (int i = 0; i < saveNodesA.n; ++i)
        saveNodesA[i] = networkA.nodes[switchIdsA[i]];

    for (int i = 0; i < saveNodesB.n; ++i)
        saveNodesB[i] = networkB.nodes[switchIdsB[i]];
    for (int i = 0; i < saveNodesT.n; ++i)
        saveNodesT[i] = networkT.nodes[switchIdsT[i]];

    // Switch and geometry optimise
    logger->debug("Switching...");
    VecF<int> optStatus_networkA(2);
    VecF<int> optStatus_SimpleGraphene(2);
    VecF<int> optStatus_TersoffGraphene(2);
    VecF<int> optStatus_TriangleRaft(2);
    VecF<int> optStatus_Bilayer(2);
    VecF<int> optStauts_BN(2);

    if (cnxType == 33) {
        // works for network version of lammps objects
        switchCnx33(switchIdsA, switchIdsB, switchIdsT);
        // works for lammps objects

        if (isSimpleGrapheneEnabled) {
            logger->debug("Switching Simple Graphene");
            SimpleGraphene.switchGraphene(switchIdsA, logger);
        }
        if (isTriangleRaftEnabled) {
            logger->debug("Switching Triangle Raft");
            Triangle_Raft.switchTriangleRaft(switchIdsA, switchIdsT, logger);
        }
        if (isBilayerEnabled)
            Bilayer.switchBilayer(switchIdsA, switchIdsT, logger);
        if (isBNEnabled) {
            logger->debug("Switching BN");
            BN.switchGraphene(switchIdsA, logger);
        }
    } else if (cnxType == 44) {
        switchCnx44(switchIdsA, switchIdsB);
    } else if (cnxType == 43) {
        switchCnx43(switchIdsA, switchIdsB);
    } else {
        logger->error("Not yet implemented (?)");
        throw std::runtime_error("Not yet implemented!");
    }
    // Rearrange nodes after switch
    bool geometryOK = true;
    geometryOK = checkThreeRingEdges(u);
    if (geometryOK)
        geometryOK = checkThreeRingEdges(v);
    if (geometryOK) {
        if (isMaintainConvexityEnabled) {
            geometryOK = convexRearrangement(cnxType, switchIdsA, switchIdsB);
            for (int i = 0; i < switchIdsA.n; ++i) {
                geometryOK = checkConvexity(switchIdsA[i]);
                if (!geometryOK)
                    break;
            }
        }
    } else {
        optStatus_SimpleGraphene = VecF<int>(3);
    }
    if (!geometryOK)
        optStatus_SimpleGraphene[0] = 4;
    int nSi;

    // Geometry optimisation of local region
    logger->debug("Optimising geometry...");
    if (geometryOK) {
        optStatus_SimpleGraphene = SimpleGraphene.globalPotentialMinimisation();
        double *localCrdsSimpleGraphene;
        double *localCrdsTersoff;
        double *localCrdsTriangleRaft;
        double *localCrdsBilayer;
        double *localCrdsBN;

        if (isSimpleGrapheneEnabled)
            localCrdsSimpleGraphene = SimpleGraphene.fetchCrds(2);
        if (isTersoffGrapheneEnabled)
            localCrdsTersoff = TersoffGraphene.fetchCrds(2);
        if (isTriangleRaftEnabled)
            localCrdsTriangleRaft = Triangle_Raft.fetchCrds(2);
        if (isBilayerEnabled)
            localCrdsBilayer = Bilayer.fetchCrds(3);
        if (isBNEnabled)
            localCrdsBN = BN.fetchCrds(2);
        if (isBilayerEnabled) {
            int natoms = Bilayer.natoms;
            int nSi = (int)(round(natoms / 3) + 0.5);
            int nO = natoms - nSi;
        }

        for (int i = 0; i < 6; ++i) {
            if (isSimpleGrapheneEnabled && isTersoffGrapheneEnabled) {
                localCrdsTersoff[2 * switchIdsA[i]] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i]] * CScaling;
                localCrdsTersoff[2 * switchIdsA[i] + 1] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * CScaling;
            }
            if (isSimpleGrapheneEnabled && isBilayerEnabled) {
                localCrdsBilayer[3 * (switchIdsA[i] * 2)] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] * 2) + 1] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] * 2 + 1)] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] * 2 + 1) + 1] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] + nSi)] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIdsA[i] + nSi) + 1] =
                    localCrdsSimpleGraphene[2 * switchIdsA[i] + 1] * SiScaling;
            }
        }
        if (isTersoffGrapheneEnabled)
            TersoffGraphene.pushCrds(2, localCrdsTersoff);
        if (isBilayerEnabled)
            Bilayer.pushCrds(3, localCrdsBilayer);

        if (isOpenMPIEnabled) {
            optStatus_networkA = localGeometryOptimisation(
                a, b, goptParamsA[1], false, isMaintainConvexityEnabled, logger);
#pragma omp parallel num_threads(3)
            {

                if (omp_get_thread_num() == 0 && isTriangleRaftEnabled) {
                    optStatus_TriangleRaft = Triangle_Raft.globalPotentialMinimisation();
                } else if (omp_get_thread_num() == 1 && isTersoffGrapheneEnabled) {
                    optStatus_TersoffGraphene =
                        TersoffGraphene.globalPotentialMinimisation();
                } else if (omp_get_thread_num() == 2 && isBilayerEnabled) {
                    optStatus_Bilayer = Bilayer.globalPotentialMinimisation();
                }
            }
        } else {

            if (isTriangleRaftEnabled) {
                optStatus_TriangleRaft = Triangle_Raft.globalPotentialMinimisation();
            }
            if (isTersoffGrapheneEnabled) {
                optStatus_TersoffGraphene =
                    TersoffGraphene.globalPotentialMinimisation();
            }
            if (isBilayerEnabled) {
                optStatus_Bilayer = Bilayer.globalPotentialMinimisation();
            }
            if (isBNEnabled) {
                optStauts_BN = BN.globalPotentialMinimisation();
            }
            optStatus_networkA = localGeometryOptimisation(
                a, b, goptParamsA[1], 0, isMaintainConvexityEnabled, logger);
        }
        if (isSimpleGrapheneEnabled)
            SimpleGrapheneEnergy = SimpleGraphene.globalPotentialEnergy();
        if (isTriangleRaftEnabled)
            TriangleRaftEnergy = Triangle_Raft.globalPotentialEnergy();
        if (isTersoffGrapheneEnabled)
            TersoffGrapheneEnergy = TersoffGraphene.globalPotentialEnergy();
        if (isBilayerEnabled)
            BilayerEnergy = Bilayer.globalPotentialEnergy();
        if (isBNEnabled)
            BNEnergy = BN.globalPotentialEnergy();
        syncCoordinates();
    } else {
        SimpleGrapheneEnergy = std::numeric_limits<double>::infinity();
        TriangleRaftEnergy = std::numeric_limits<double>::infinity();
        TersoffGrapheneEnergy = std::numeric_limits<double>::infinity();
        BilayerEnergy = std::numeric_limits<double>::infinity();
        BNEnergy = std::numeric_limits<double>::infinity();
    }
    logger->info("Accepting or rejecting...");
    bool isAccepted = false;
    if (mcRoutine == 1)
        isAccepted = mc.acceptanceCriterion(SimpleGrapheneEnergy,
                                            saveEnergySimpleGraphene, 1.0);
    else if (mcRoutine == 2)
        isAccepted = mc.acceptanceCriterion(TriangleRaftEnergy,
                                            saveEnergyTriangleRaft, 7.3448);
    else if (mcRoutine == 5)
        isAccepted = mc.acceptanceCriterion(BNEnergy, saveEnergyBN, 7.0);
    if (isAccepted) {
        if (mcRoutine == 1)
            logger->info("Accepted Move, Ei = {} Ef = {}", saveEnergySimpleGraphene,
                         SimpleGrapheneEnergy);
        else if (mcRoutine == 2)
            logger->info("Accepted Move, Ei = {} Ef = {}", saveEnergyTriangleRaft,
                         TriangleRaftEnergy);
        else if (mcRoutine == 5)
            logger->info("Accepted Move, Ei = {} Ef = {}", saveEnergyBN, BNEnergy);
        else
            logger->info("MC ROUTINE : {}", mcRoutine);
        logger->info("Syncing coordinates...");
        syncCoordinates();
    } else {
        if (mcRoutine == 1)
            logger->info("Rejected Move, Ei = {} Ef = {}", saveEnergySimpleGraphene,
                         SimpleGrapheneEnergy);
        else if (mcRoutine == 2)
            logger->info("Rejected Move, Ei = {} Ef = {}", saveEnergyTriangleRaft,
                         TriangleRaftEnergy);
        else if (mcRoutine == 5)
            logger->info("Rejected Move, Ei = {} Ef = {}", saveEnergyBN, BNEnergy);

        crds = saveCrds;
        networkA.nodeDistribution = saveNodeDistA;
        networkA.edgeDistribution = saveEdgeDistA;
        networkB.nodeDistribution = saveNodeDistB;
        networkB.edgeDistribution = saveEdgeDistB;

        for (int i = 0; i < saveNodesA.n; ++i)
            networkA.nodes[switchIdsA[i]] = saveNodesA[i];
        for (int i = 0; i < saveNodesB.n; ++i)
            networkB.nodes[switchIdsB[i]] = saveNodesB[i];
        for (int i = 0; i < saveNodesT.n; ++i)
            networkT.nodes[switchIdsT[i]] = saveNodesT[i];

        if (isSimpleGrapheneEnabled)
            SimpleGrapheneEnergy = saveEnergySimpleGraphene;
        if (isTriangleRaftEnabled)
            TriangleRaftEnergy = saveEnergyTriangleRaft;
        if (isTersoffGrapheneEnabled)
            TersoffGrapheneEnergy = saveEnergyTersoffGraphene;
        if (isBilayerEnabled)
            BilayerEnergy = saveEnergyBilayer;
        if (isBNEnabled)
            BNEnergy = saveEnergyBN;
        if (isSimpleGrapheneEnabled)
            SimpleGraphene.pushCrds(2, saveCrdsSimpleGraphene);
        if (isTriangleRaftEnabled)
            Triangle_Raft.pushCrds(2, saveCrdsTriangleRaft);
        if (isBilayerEnabled)
            Bilayer.pushCrds(3, saveCrdsBilayer);
        if (isTersoffGrapheneEnabled)
            TersoffGraphene.pushCrds(2, saveCrdsTersoffGraphene);
        if (isBNEnabled)
            BN.pushCrds(2, saveCrdsBN);

        if (isSimpleGrapheneEnabled) {
            SimpleGraphene.revertGraphene(switchIdsA, logger);
            SimpleGraphene.globalPotentialMinimisation();
        }
        if (isTriangleRaftEnabled) {
            Triangle_Raft.revertTriangleRaft(switchIdsA, switchIdsT, logger);
            Triangle_Raft.globalPotentialMinimisation();
        }
        if (isBilayerEnabled) {
            Bilayer.revertBilayer(switchIdsA, switchIdsT, logger);
            Bilayer.globalPotentialMinimisation();
        }
        if (isBNEnabled) {
            BN.revertGraphene(switchIdsA, logger);
            BN.globalPotentialMinimisation();
        }
    }
    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force)
     * 2=unsuccessful(it limit) 3=unsuccessful(intersection)
     * 4=unsuccessful(non-convex) [2] optimisation iterations */
    VecF<int> status(3);
    status[0] = isAccepted;
    status[1] = optStatus_SimpleGraphene[0];
    status[2] = optStatus_SimpleGraphene[1];
    return status;
}

// Single monte carlo switching move
VecF<int> LinkedNetwork::monteCarloSwitchMove(Network network, double &energy,
                                              LoggerPtr logger) {

    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */

    // Select valid random connection - that will not violate connection limits
    int a;
    int b;
    int u;
    int v;
    VecF<int> switchIdsA;
    VecF<int> switchIdsB;
    VecF<int> switchIdsT;
    bool foundValidMove = false;
    int cnxType;
    for (int i = 0; i < networkA.nodes.n * networkA.nodes.n;
         ++i) { // catch in case cannot find any valid moves
        cnxType = pickRandomConnection(a, b, u, v, mtGen, SelectionType::RANDOM, logger);
        foundValidMove = generateSwitchIDs(cnxType, switchIdsA, switchIdsB,
                                           switchIdsT, a, b, u, v, logger);
        if (foundValidMove)
            break;
    }
    if (!foundValidMove) {
        logger->critical("Cannot find any valid switch moves");
        throw std::runtime_error("Cannot find any valid switch moves");
    }

    // Save current state
    double saveEnergy = energy;
    VecF<double> saveCrds = crds;
    VecF<int> saveNodeDistA = networkA.nodeDistribution;
    VecF<int> saveNodeDistB = networkB.nodeDistribution;
    VecF<VecF<int>> saveEdgeDistA = networkA.edgeDistribution;
    VecF<VecF<int>> saveEdgeDistB = networkB.edgeDistribution;
    VecF<Node> saveNodesA(switchIdsA.n);
    VecF<Node> saveNodesB(switchIdsB.n);
    VecF<Node> saveNodesT(switchIdsT.n);
    for (int i = 0; i < saveNodesA.n; ++i)
        saveNodesA[i] = networkA.nodes[switchIdsA[i]];
    for (int i = 0; i < saveNodesB.n; ++i)
        saveNodesB[i] = networkB.nodes[switchIdsB[i]];
    for (int i = 0; i < saveNodesT.n; ++i)
        saveNodesT[i] = networkT.nodes[switchIdsT[i]];

    // Switch and geometry optimise
    logger->info("Switching...");
    VecF<int> optStatus(3);
    if (cnxType == 33)
        switchCnx33(switchIdsA, switchIdsB, switchIdsT);
    else if (cnxType == 44)
        switchCnx44(switchIdsA, switchIdsB);
    else if (cnxType == 43)
        switchCnx43(switchIdsA, switchIdsB);
    else {
        logger->critical("Not yet implemented cnxType: {}", cnxType);
        throw std::runtime_error("Not yet implemented");
    }
    // Rearrange nodes after switch
    bool geometryOK = true;
    geometryOK = checkThreeRingEdges(u);
    if (geometryOK)
        geometryOK = checkThreeRingEdges(v);
    if (geometryOK) {
        if (isMaintainConvexityEnabled) {
            geometryOK = convexRearrangement(cnxType, switchIdsA, switchIdsB);
            for (int i = 0; i < switchIdsA.n; ++i) {
                geometryOK = checkConvexity(switchIdsA[i]);
                if (!geometryOK)
                    break;
            }
        }
    } else {
        optStatus = VecF<int>(3);
    }
    if (!geometryOK)
        optStatus[0] = 4;

    // Geometry optimisation of local region
    logger->info("Optimising geometry...");
    if (geometryOK) {

        optStatus = globalGeometryOptimisation(false, false, network, logger);
        energy =
            globalPotentialEnergy(0, isMaintainConvexityEnabled, network, logger);
    } else
        energy = std::numeric_limits<double>::infinity();

    // Accept or reject
    bool isAccepted = mc.acceptanceCriterion(energy, saveEnergy, 1.00);
    if (isAccepted) {
        logger->info("Accepted MC Move Ei = {} Ef = {}", saveEnergy, energy);
        energy = saveEnergy;
        crds = saveCrds;
        networkA.nodeDistribution = saveNodeDistA;
        networkA.edgeDistribution = saveEdgeDistA;
        networkB.nodeDistribution = saveNodeDistB;
        networkB.edgeDistribution = saveEdgeDistB;
        for (int i = 0; i < saveNodesA.n; ++i)
            networkA.nodes[switchIdsA[i]] = saveNodesA[i];
        for (int i = 0; i < saveNodesB.n; ++i)
            networkB.nodes[switchIdsB[i]] = saveNodesB[i];
        for (int i = 0; i < saveNodesT.n; ++i)
            networkT.nodes[switchIdsT[i]] = saveNodesT[i];
    } else {
        logger->info("Rejected MC Move Ei = {} Ef = {}", saveEnergy, energy);
        syncCoordinates();
    }

    /* Status report
     * [0] accepted/rejected 1/0
     * [1] optimisation code 0=successful 1=successful(zero force)
     * 2=unsuccessful(it limit) 3=unsuccessful(intersection)
     * 4=unsuccessful(non-convex) [2] optimisation iterations */
    VecF<int> status(3);
    status[0] = isAccepted;
    status[1] = optStatus[0];
    status[2] = optStatus[1];
    return status;
}

// Calculate potential energy of entire system
double LinkedNetwork::globalPotentialEnergy(bool useIntx, bool isMaintainConvexityEnabledArg,
                                            Network network, LoggerPtr logger) {

    /* Potential model
     * Bonds as harmonic
     * Angles as harmonic
     * Local line intersections */

    // Set up potential model
    // Bond and angle harmonics
    VecR<int> bonds(0, network.nodes.n * 6);
    VecR<int> angles(0, network.nodes.n * network.maxNetCnxs * 3);
    VecR<double> bondParams(0, network.nodes.n * 6);
    VecR<double> angleParams(0, network.nodes.n * network.maxNetCnxs * 3);
    if (network.nodes.n > networkA.nodes.n) {
        for (int i = 0; i < network.nodes.n; ++i) {
            generateHarmonicsOnly(i, bonds, bondParams, network);
        }
    } else {
        for (int i = 0; i < network.nodes.n; ++i) {
            generateHarmonics(i, bonds, bondParams, angles, angleParams, network,
                              logger);
        }
    }
    // Intersections
    VecR<int> intersections(0, network.nodes.n * 1000);
    if (useIntx) {
        for (int i = 0; i < networkB.nodes.n; ++i) {
            generateRingIntersections(i, intersections);
        }
    }

    // Assign to fixed size arrays
    VecF<int> bnds(bonds.n);
    VecF<int> angs(angles.n);
    VecF<int> intx(intersections.n);
    VecF<double> bndP(bondParams.n);
    VecF<double> angP(angleParams.n);
    VecF<double> intxP; // placeholder for intersections
    for (int i = 0; i < bnds.n; ++i)
        bnds[i] = bonds[i];
    for (int i = 0; i < bndP.n; ++i)
        bndP[i] = bondParams[i];
    for (int i = 0; i < angs.n; ++i)
        angs[i] = angles[i];
    for (int i = 0; i < angP.n; ++i)
        angP[i] = angleParams[i];
    for (int i = 0; i < intersections.n; ++i)
        intx[i] = intersections[i];

    // Potential model based on geometry code
    double potEnergy = 0.0;
    if (network.geometryCode == "2DE") {
        if (!isMaintainConvexityEnabledArg) {
            HI2DP potModel(network.dimensions[0], network.dimensions[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx)
                potModel.setIntersections(intx, intxP);
            potEnergy = potModel.function(crds);
        } else {
            HRI2DP potModel(network.dimensions[0], network.dimensions[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx)
                potModel.setIntersections(intx, intxP);
            potEnergy = potModel.function(crds);
        }
    } else if (network.geometryCode == "2DS") {
        VecF<int> constrained(network.nodes.n);
        for (int i = 0; i < network.nodes.n; ++i)
            constrained[i] = i;
        HI3DS potModel;
        potModel.setBonds(bnds, bndP);
        potModel.setAngles(angs, angP);
        potModel.setGeomConstraints(constrained, potParamsC);
        if (useIntx)
            potModel.setIntersections(intx, intxP);
        potEnergy = potModel.function(crds);
    }

    else if (network.geometryCode == "2DEtr") {
        HLJ2DP potModel(network.dimensions[0], network.dimensions[1]);
        potModel.setBonds(bnds, bndP);
        potModel.setAngles(angs, angP);
        if (useIntx)
            potModel.setIntersections(intx, intxP);
        potEnergy = potModel.function(crds);
    } else {
        logger->critical("Not yet implemented geometry code: {}",
                         network.geometryCode);
        throw std::runtime_error("Not yet implemented");
    }
    // Convexity
    if (isMaintainConvexityEnabledArg) {
        bool convex = checkConvexity();
        if (!convex)
            potEnergy = std::numeric_limits<double>::infinity();
    }
    return potEnergy;
}

// Geometry optimise entire system
VecF<int>
LinkedNetwork::globalGeometryOptimisation(bool useIntx,
                                          bool isMaintainConvexityEnabledArg,
                                          Network network, LoggerPtr logger) {
    auto start_GEO = std::chrono::system_clock::now();
    /* Potential model
     * Bonds as harmonic
     * Angles as approximated harmonic
     * Local line intersections */
    int maxCnxs = network.maxNetCnxs;
    // Set up potential model
    // Bond and angle harmonics
    VecR<int> bonds(0, network.nodes.n * maxCnxs * 2);
    VecR<int> angles(0, network.nodes.n * maxCnxs * 3);
    VecR<double> bondParams(0, network.nodes.n * maxCnxs * 3);
    VecR<double> angleParams(0, network.nodes.n * maxCnxs * 3);
    if (network.nodes.n > networkA.nodes.n) {
        for (int i = 0; i < network.nodes.n; ++i) {
            generateHarmonicsOnly(i, bonds, bondParams, network);
        }
    } else {
        for (int i = 0; i < network.nodes.n; ++i) {
            generateHarmonics(i, bonds, bondParams, angles, angleParams, network,
                              logger);
        }
    }
    // Intersections
    VecR<int> intersections(0, network.nodes.n * 1000);

    // Assign to fixed size arrays
    VecF<int> bnds(bonds.n);
    VecF<int> angs(angles.n);
    VecF<int> intx(intersections.n);
    VecF<double> bndP(bondParams.n);
    VecF<double> angP(angleParams.n);
    VecF<double> intxP; // placeholder for intersections
    for (int i = 0; i < bnds.n; ++i)
        bnds[i] = bonds[i];
    for (int i = 0; i < bndP.n; ++i)
        bndP[i] = bondParams[i];
    for (int i = 0; i < angs.n; ++i)
        angs[i] = angles[i];
    for (int i = 0; i < angP.n; ++i)
        angP[i] = angleParams[i];
    for (int i = 0; i < intersections.n; ++i)
        intx[i] = intersections[i];

    VecF<int> optStatus(2);

    // Geometry optimise
    if (network.geometryCode == "2DE") {
        if (!isMaintainConvexityEnabledArg) {
            HI2DP potModel(network.dimensions[0], network.dimensions[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx)
                potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) <
                std::numeric_limits<double>::infinity()) { // only optimise if no line
                                                           // intersections
                SteepestDescentArmijoMultiDim<HI2DP> optimiser(
                    goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;
            }
        } else {
            HRI2DP potModel(network.dimensions[0], network.dimensions[1]);
            potModel.setBonds(bnds, bndP);
            potModel.setAngles(angs, angP);
            if (useIntx)
                potModel.setIntersections(intx, intxP);
            if (potModel.function(crds) <
                std::numeric_limits<double>::infinity()) { // only optimise if no line
                                                           // intersections
                SteepestDescentArmijoMultiDim<HRI2DP> optimiser(goptParamsA[0], goptParamsB[0], goptParamsB[1]);
                optStatus = optimiser(potModel, crds);
            } else {
                optStatus[0] = 3;
                optStatus[1] = 0;
            }
        }
    } else if (network.geometryCode == "2DEtr") {
        logger->info("Minimising Triangle Raft...");
        HLJ2DP potModel(network.dimensions[0], network.dimensions[1]);
        potModel.setBonds(bnds, bndP);
        if (potModel.function(crds) <
            std::numeric_limits<double>::infinity()) { // only optimise if no line
                                                       // intersections
            SteepestDescentArmijoMultiDim<HLJ2DP> optimiser(
                goptParamsA[0], goptParamsB[0], goptParamsB[1]);
            optStatus = optimiser(potModel, crds);
        } else {
            optStatus[0] = 3;
            optStatus[1] = 0;
        }
    }
    auto end_GEO = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_GEO = end_GEO - start_GEO;
    logger->info("Global geometry optimisation took {} seconds",
                 elapsed_GEO.count());
    return optStatus;
}

// Geometry optimise subsection of system by only including interactions in a
// specified range
VecF<int> LinkedNetwork::localGeometryOptimisation(
    int centreA, int centreB, int extent, bool useIntx,
    bool isMaintainConvexityEnabledArg, LoggerPtr logger) {
    /* Find three regions
     * 1) local (full interactions)
     * 2) fixed inner (interactions with local and fixed but immobile)
     * 3) fixed outer (interactions with fixed inner but immobile) */
    VecR<int> local;
    VecR<int> fixedInner;
    VecR<int> fixedOuter;
    networkA.findLocalRegion(centreA, centreB, extent, local, fixedInner,
                             fixedOuter);

    // Harmonics
    int reserveSize = (local.n + fixedInner.n) * maxACnxs * 3;
    VecR<int> bonds(0, reserveSize);
    VecR<int> angles(0, reserveSize);
    VecR<double> bondParams(0, reserveSize);
    VecR<double> angleParams(0, reserveSize);
    for (int i = 0; i < local.n; ++i) {
        generateHarmonics(local[i], bonds, bondParams, angles, angleParams,
                          networkA, logger);
    }
    for (int i = 0; i < fixedInner.n; ++i) {
        generateHarmonics(fixedInner[i], bonds, bondParams, angles, angleParams,
                          networkA, logger);
    }

    // Intersections - Expensive, turn off and use in global potential energy
    VecR<int> intersections(0, local.n * 1000);

    // Assign to fixed size arrays
    VecF<int> bnds(bonds.n);
    VecF<int> fixd(fixedInner.n + fixedOuter.n);
    VecF<int> angs(angles.n);
    VecF<int> intx(intersections.n);
    VecF<double> bndP(bondParams.n);
    VecF<double> angP(angleParams.n);
    VecF<double> intxP; // placeholder for intersections
    for (int i = 0; i < bnds.n; ++i)
        bnds[i] = bonds[i];
    for (int i = 0; i < bndP.n; ++i)
        bndP[i] = bondParams[i];
    for (int i = 0; i < angs.n; ++i)
        angs[i] = angles[i];
    for (int i = 0; i < angP.n; ++i)
        angP[i] = angleParams[i];
    for (int i = 0; i < fixedInner.n; ++i)
        fixd[i] = fixedInner[i];
    for (int i = 0; i < fixedOuter.n; ++i)
        fixd[i + fixedInner.n] = fixedOuter[i];
    for (int i = 0; i < intersections.n; ++i)
        intx[i] = intersections[i];

    // Geometry optimise
    VecF<int> optStatus(2);
    if (!isMaintainConvexityEnabledArg) {
        HI2DP potModel(networkA.dimensions[0], networkA.dimensions[1]);
        potModel.setBonds(bnds, bndP);
        potModel.setAngles(angs, angP);
        potModel.setFixedAtoms(fixd);
        if (useIntx)
            potModel.setIntersections(intx, intxP);
        if (potModel.function(crds) <
            std::numeric_limits<double>::infinity()) { // optimise if no line
                                                       // intersections
            SteepestDescentArmijoMultiDim<HI2DP> optimiser(
                goptParamsA[0], goptParamsB[0], goptParamsB[1]);
            optStatus = optimiser(potModel, crds);
        } else {
            optStatus[0] = 3;
            optStatus[1] = 0;
        }
    } else {
        HRI2DP potModel(networkA.dimensions[0], networkA.dimensions[1]);
        potModel.setBonds(bnds, bndP);
        potModel.setAngles(angs, angP);
        potModel.setFixedAtoms(fixd);
        if (useIntx)
            potModel.setIntersections(intx, intxP);
        if (potModel.function(crds) <
            std::numeric_limits<double>::infinity()) { // optimise if no line
                                                       // intersections
            SteepestDescentArmijoMultiDim<HRI2DP> optimiser(
                goptParamsA[0], goptParamsB[0], goptParamsB[1]);
            optStatus = optimiser(potModel, crds);
        } else {
            optStatus[0] = 3;
            optStatus[1] = 0;
        }
    }
    return optStatus;
}

// Generate harmonic potentials corresponding to bonds and angles associated
// with a given node in lattice A
void LinkedNetwork::generateHarmonics(int id, VecR<int> &bonds,
                                      VecR<double> &bondParams,
                                      VecR<int> &angles,
                                      VecR<double> &angleParams,
                                      Network network, LoggerPtr logger) {
    // Potential parameters
    double bk = potParamsB[0];
    double br = potParamsB[1]; // bond k and r0
    double br0 = sqrt(3) * br;
    double ak;
    double act; // angle k, cos theta0
    // Harmonics
    int cnd; // coordination
    int id0 = id;
    int id1;
    int id2;
    cnd = network.nodes[id0].netCnxs.n;
    ak = potParamsA[0];
    act = cos(2.0 * M_PI / cnd);
    for (int i = 0; i < cnd; ++i) {
        id1 = network.nodes[id0].netCnxs[i];
        id2 = network.nodes[id0].netCnxs[(i + 1) % cnd];
        // bonds
        if (id0 < id1) {
            // prevent double counting
            if (id0 >= networkA.nodes.n && id1 >= networkA.nodes.n) {
                bonds.addValue(id0);
                bonds.addValue(id1);
                bondParams.addValue(bk);
                bondParams.addValue(br0);
            } else {
                bonds.addValue(id0);
                bonds.addValue(id1);
                bondParams.addValue(bk);
                bondParams.addValue(br);
            }
        }
        // angles
        angles.addValue(id0);
        angles.addValue(id1);
        angles.addValue(id2);
        angleParams.addValue(ak);
        angleParams.addValue(act);
    }
}

// Generate harmonic potentials corresponding to bonds and angles associated
// with a given node in lattice A
void LinkedNetwork::generateHarmonicsOnly(int id, VecR<int> &bonds,
                                          VecR<double> &bondParams,
                                          Network network) {

    // Potential parameters
    double bk = potParamsB[0];
    double br = potParamsB[1]; // bond k and r0
    double br0 = sqrt(3) * br;

    // Harmonics
    int cnd; // coordination
    int id0 = id;
    int id1;
    int id2;
    cnd = network.nodes[id0].netCnxs.n;

    int sio_bond_count = 0;
    int o_o_bond_count = 0;

    for (int i = 0; i < cnd; ++i) {
        id1 = network.nodes[id0].netCnxs[i];
        id2 = network.nodes[id0].netCnxs[(i + 1) % cnd];
        // bonds
        if (id0 < id1) {
            if (network.nodes.n > networkA.nodes.n) {
                // prevent double counting
                if (id0 >= networkA.nodes.n && id1 >= networkA.nodes.n) {
                    bonds.addValue(id0);
                    bonds.addValue(id1);
                    bondParams.addValue(bk);
                    bondParams.addValue(br0 / 2);
                    sio_bond_count += 1;
                } else {
                    bonds.addValue(id0);
                    bonds.addValue(id1);
                    bondParams.addValue(bk);
                    bondParams.addValue(br / 2);
                    o_o_bond_count += 1;
                }
            } else {
                bonds.addValue(id0);
                bonds.addValue(id1);
                bondParams.addValue(bk);
                bondParams.addValue(br);
            }
        }
    }
}

// Generate ring edge intersections for a specific ring
void LinkedNetwork::generateRingIntersections(int rId,
                                              VecR<int> &intersections) {

    int rCnd = networkB.nodes[rId].netCnxs.n;
    int nCnd = networkB.nodes[rId].dualCnxs.n;
    int e0;
    int e1;
    int e2;
    int e3;
    for (int i = 0; i < rCnd; ++i) { // loop over neighbouring rings
        int rId0 = networkB.nodes[rId].netCnxs[i];
        if (rId < rId0) {                    // prevent double counting
            for (int j = 0; j < nCnd; ++j) { // loop over nodes
                e0 = networkB.nodes[rId].dualCnxs[j];
                e1 = networkB.nodes[rId].dualCnxs[(j + 1) % nCnd];
                int nCnd0 = networkB.nodes[rId0].dualCnxs.n;
                for (int k = 0; k < nCnd0; ++k) {
                    e2 = networkB.nodes[rId0].dualCnxs[k];
                    e3 = networkB.nodes[rId0].dualCnxs[(k + 1) % nCnd0];
                    if (e0 != e2 && e0 != e3 && e1 != e2 && e1 != e3) {
                        intersections.addValue(e0);
                        intersections.addValue(e1);
                        intersections.addValue(e2);
                        intersections.addValue(e3);
                    }
                }
            }
        }
    }
}

// Generate intersections required to maintain convexity for a given node
void LinkedNetwork::generateConvexIntersections(int nId,
                                                VecR<int> &intersections) {

    int cnd = networkA.nodes[nId].netCnxs.n;
    int id0;
    int id1;
    int id2;
    for (int i = 0; i < cnd; ++i) {
        id0 = networkA.nodes[nId].netCnxs[i];
        id1 = networkA.nodes[nId].netCnxs[(i + 1) % cnd];
        id2 = networkA.nodes[nId].netCnxs[(i + 2) % cnd];
        intersections.addValue(id0);
        intersections.addValue(id1);
        intersections.addValue(nId);
        intersections.addValue(id2);
    }
}

// Update networks with geometry optimised coordinates
void LinkedNetwork::syncCoordinates() {

    // Sync T coordinates
    for (int i = 0; i < networkA.nodes.n; ++i) {
        networkT.nodes[i].crd[0] = crds[2 * i];
        networkT.nodes[i].crd[1] = crds[2 * i + 1];
    }
    // Sync A coordinates
    if (networkA.geometryCode == "2DE") {
        for (int i = 0; i < networkA.nodes.n; ++i) {
            networkA.nodes[i].crd[0] = crds[2 * i];
            networkA.nodes[i].crd[1] = crds[2 * i + 1];
        }
    } else {
        for (int i = 0; i < networkA.nodes.n; ++i) {
            networkA.nodes[i].crd[0] = crds[3 * i];
            networkA.nodes[i].crd[1] = crds[3 * i + 1];
            networkA.nodes[i].crd[2] = crds[3 * i + 2];
        }
    }

    // Sync B coordinates
    for (int i = 0; i < networkB.nodes.n; ++i) {
        VecF<double> x(networkB.nodes[i].dualCnxs.n);
        VecF<double> y(networkB.nodes[i].dualCnxs.n);
        for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
            x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
            y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
        }
        VecF<double> origin(2);
        origin[0] = x[0];
        origin[1] = y[0];
        x -= origin[0];
        y -= origin[1];
        for (int j = 0; j < x.n; ++j)
            x[j] -= networkB.dimensions[0] * nearbyint(x[j] * networkB.rpb[0]);
        for (int j = 0; j < y.n; ++j)
            y[j] -= networkB.dimensions[1] * nearbyint(y[j] * networkB.rpb[1]);
        VecF<double> c(2);
        c[0] = origin[0] + vMean(x);
        c[1] = origin[1] + vMean(y);
        networkB.nodes[i].crd = c;
    }
}

// Update networks with geometry optimised coordinates
void LinkedNetwork::syncCoordinatesTD() {

    // Sync A coordinates
    for (int i = 0; i < networkA.nodes.n; ++i) {
        networkA.nodes[i].crd[0] = crds[2 * i];
        networkA.nodes[i].crd[1] = crds[2 * i + 1];
    }
    for (int i = 0; i < networkT.nodes.n; ++i) {
        networkT.nodes[i].crd[0] = crds[2 * i];
        networkT.nodes[i].crd[1] = crds[2 * i + 1];
    }

    // Sync B coordinates
    for (int i = 0; i < networkB.nodes.n; ++i) {
        VecF<double> x(networkB.nodes[i].dualCnxs.n);
        VecF<double> y(networkB.nodes[i].dualCnxs.n);
        for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
            x[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
            y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
        }
        VecF<double> origin(2);
        origin[0] = x[0];
        origin[1] = y[0];
        x -= origin[0];
        y -= origin[1];
        for (int j = 0; j < x.n; ++j)
            x[j] -= networkB.dimensions[0] * nearbyint(x[j] * networkB.rpb[0]);
        for (int j = 0; j < y.n; ++j)
            y[j] -= networkB.dimensions[1] * nearbyint(y[j] * networkB.rpb[1]);
        VecF<double> c(2);
        c[0] = origin[0] + vMean(x);
        c[1] = origin[1] + vMean(y);
        networkB.nodes[i].crd = c;
    }
}

// Get normalised probability distribution of nodes of each size in given
// lattice
VecF<double> LinkedNetwork::getNodeDistribution(std::string_view lattice) {
    if (lattice == "A")
        return networkA.getNodeDistribution();
    else
        return networkB.getNodeDistribution();
}

// Get unnormalised probability distribution of node connections in given
// lattice
VecF<VecF<int>> LinkedNetwork::getEdgeDistribution(std::string_view lattice) {

    if (lattice == "A")
        return networkA.edgeDistribution;
    else
        return networkB.edgeDistribution;
}

// Get Aboav-Weaire fitting parameters
VecF<double> LinkedNetwork::getAboavWeaire(std::string_view lattice) {

    if (lattice == "A")
        return networkA.aboavWeaireParams();
    else
        return networkB.aboavWeaireParams();
}

// Get assortativity
double LinkedNetwork::getAssortativity(std::string_view lattice) {

    if (lattice == "A")
        return networkA.assortativity();
    else
        return networkB.assortativity();
}

// Get alpha estimate
double LinkedNetwork::getAboavWeaireEstimate(std::string_view lattice) {

    if (lattice == "A")
        return networkA.aboavWeaireEstimate();
    else
        return networkB.aboavWeaireEstimate();
}

// Get entropy
VecF<double> LinkedNetwork::getEntropy(std::string_view lattice) {

    if (lattice == "A")
        return networkA.entropy();
    else
        return networkB.entropy();
}

// Get cluster information
double LinkedNetwork::getMaxCluster(std::string_view lattice, int nodeCnd) {

    if (lattice == "A")
        return networkA.maxCluster(nodeCnd);
    else
        return networkB.maxCluster(nodeCnd);
}

// Get cluster information
VecF<int> LinkedNetwork::getMaxClusters(std::string_view lattice, int minCnd,
                                        int maxCnd) {

    if (lattice == "A")
        return networkA.maxClusters(minCnd, maxCnd, 3, 2);
    else
        return networkB.maxClusters(minCnd, maxCnd, 3, 2);
}

// Check linked networks for consistency
bool LinkedNetwork::checkConsistency() {

    bool checkCnx = checkCnxConsistency();
    bool checkDesc = checkDescriptorConsistency();
    bool check = checkCnx * checkDesc;

    return check;
}

// Check linked networks have mutual network and dual connections
bool LinkedNetwork::checkCnxConsistency() {

    // check number of network connections is equal to number of dual connections
    bool netDualEquality = true;
    for (int i = 0; i < networkA.nodes.n; ++i) {
        if (networkA.nodes[i].netCnxs.n != networkA.nodes[i].dualCnxs.n)
            netDualEquality = false;
    }
    for (int i = 0; i < networkB.nodes.n; ++i) {
        if (networkB.nodes[i].netCnxs.n != networkB.nodes[i].dualCnxs.n)
            netDualEquality = false;
    }

    // check mutual network connections
    bool mutualNetCnx = true;
    int id0;
    int id1;
    for (int i = 0; i < networkA.nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < networkA.nodes[i].netCnxs.n; ++j) {
            id1 = networkA.nodes[i].netCnxs[j];
            mutualNetCnx = vContains(networkA.nodes[id1].netCnxs, id0);
        }
    }
    for (int i = 0; i < networkB.nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < networkB.nodes[i].netCnxs.n; ++j) {
            id1 = networkB.nodes[i].netCnxs[j];
            mutualNetCnx = vContains(networkB.nodes[id1].netCnxs, id0);
        }
    }

    // check mutual dual connections
    bool mutualDualCnx = true;
    for (int i = 0; i < networkA.nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < networkA.nodes[i].dualCnxs.n; ++j) {
            id1 = networkA.nodes[i].dualCnxs[j];
            mutualDualCnx = vContains(networkB.nodes[id1].dualCnxs, id0);
        }
    }
    for (int i = 0; i < networkB.nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
            id1 = networkB.nodes[i].dualCnxs[j];
            mutualDualCnx = vContains(networkA.nodes[id1].dualCnxs, id0);
        }
    }

    // check network connections are neighbours by lying on same ring (some highly
    // strained cases could give a false positive)
    bool nbNetCnx = true;
    for (int i = 0; i < networkA.nodes.n; ++i) {
        int nCnxs = networkA.nodes[i].netCnxs.n;
        for (int j = 0; j < nCnxs; ++j) {
            id0 = networkA.nodes[i].netCnxs[j];
            id1 = networkA.nodes[i].netCnxs[(j + 1) % nCnxs];
            VecR<int> common = vCommonValues(networkA.nodes[id0].dualCnxs,
                                             networkA.nodes[id1].dualCnxs);
            if (common.n == 0)
                nbNetCnx = false;
        }
    }
    for (int i = 0; i < networkB.nodes.n; ++i) {
        int nCnxs = networkB.nodes[i].netCnxs.n;
        for (int j = 0; j < nCnxs; ++j) {
            id0 = networkB.nodes[i].netCnxs[j];
            id1 = networkB.nodes[i].netCnxs[(j + 1) % nCnxs];
            VecR<int> common = vCommonValues(networkB.nodes[id0].dualCnxs,
                                             networkB.nodes[id1].dualCnxs);
            if (common.n == 0)
                nbNetCnx = false;
        }
    }

    // check dual connections are neighbours by lying on same ring (some highly
    // strained cases could give a false positive)
    bool nbDualCnx = true;
    for (int i = 0; i < networkA.nodes.n; ++i) {
        int nCnxs = networkA.nodes[i].dualCnxs.n;
        for (int j = 0; j < nCnxs; ++j) {
            id0 = networkA.nodes[i].dualCnxs[j];
            id1 = networkA.nodes[i].dualCnxs[(j + 1) % nCnxs];
            VecR<int> common = vCommonValues(networkB.nodes[id0].dualCnxs,
                                             networkB.nodes[id1].dualCnxs);
            common.delValue(i);
            if (common.n == 0)
                nbDualCnx = false;
        }
    }
    for (int i = 0; i < networkB.nodes.n; ++i) {
        int nCnxs = networkB.nodes[i].dualCnxs.n;
        for (int j = 0; j < nCnxs; ++j) {
            id0 = networkB.nodes[i].dualCnxs[j];
            id1 = networkB.nodes[i].dualCnxs[(j + 1) % nCnxs];
            VecR<int> common = vCommonValues(networkA.nodes[id0].dualCnxs,
                                             networkA.nodes[id1].dualCnxs);
            common.delValue(i);
            if (common.n == 0)
                nbDualCnx = false;
        }
    }

    // check expected number of auxiliary connections
    bool numAux = true;
    for (int i = 0; i < networkB.nodes.n; ++i) {
        int expAux = 0;
        for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j)
            expAux += networkA.nodes[networkB.nodes[i].dualCnxs[j]].netCnxs.n - 3;
        if (networkB.nodes[i].auxCnxs.n != expAux)
            numAux = false;
    }
    numAux = true;

    // check mutual auxiliary connections
    bool mutualAuxCnx = true;
    for (int i = 0; i < networkB.nodes.n; ++i) {
        id0 = i;
        for (int j = 0; j < networkB.nodes[i].auxCnxs.n; ++j) {
            id1 = networkB.nodes[i].auxCnxs[j];
            mutualAuxCnx = vContains(networkB.nodes[id1].auxCnxs, id0);
        }
    }
    mutualAuxCnx = true;

    // overall flag
    bool consistent = netDualEquality * mutualNetCnx * mutualDualCnx * nbNetCnx *
                      nbDualCnx * numAux * mutualAuxCnx;

    return consistent;
}

// Check linked networks have accurate descriptors
bool LinkedNetwork::checkDescriptorConsistency() {

    VecF<int> checkNodeA(networkA.nodeDistribution);
    VecF<int> checkNodeB(networkB.nodeDistribution);
    VecF<VecF<int>> checkEdgeA(networkA.edgeDistribution);
    VecF<VecF<int>> checkEdgeB(networkB.edgeDistribution);

    // Check node distribution
    bool nodeA;
    bool nodeB;
    checkNodeA = 0;
    checkNodeB = 0;
    for (int i = 0; i < networkA.nodes.n; ++i) {
        int n = networkA.nodes[i].netCnxs.n;
        ++checkNodeA[n];
    }
    for (int i = 0; i < networkB.nodes.n; ++i) {
        int n = networkB.nodes[i].netCnxs.n;
        ++checkNodeB[n];
    }
    nodeA = checkNodeA == networkA.nodeDistribution;
    nodeB = checkNodeB == networkB.nodeDistribution;

    // Check edge distribution
    bool edgeA = true;
    bool edgeB = true;
    for (int i = 0; i < checkEdgeA.n; ++i)
        checkEdgeA[i] = 0;
    for (int i = 0; i < checkEdgeB.n; ++i)
        checkEdgeB[i] = 0;
    for (int i = 0; i < networkA.nodes.n; ++i) {
        int m = networkA.nodes[i].netCnxs.n;
        for (int j = 0; j < m; ++j) {
            int n = networkA.nodes[networkA.nodes[i].netCnxs[j]].netCnxs.n;
            ++checkEdgeA[m][n];
        }
    }
    for (int i = 0; i < networkB.nodes.n; ++i) {
        int m = networkB.nodes[i].netCnxs.n;
        for (int j = 0; j < m; ++j) {
            int n = networkB.nodes[networkB.nodes[i].netCnxs[j]].netCnxs.n;
            ++checkEdgeB[m][n];
        }
    }
    for (int i = 0; i < checkEdgeA.n; ++i) {
        if (!(checkEdgeA[i] == networkA.edgeDistribution[i])) {
            edgeA = false;
            break;
        }
    };
    for (int i = 0; i < checkEdgeB.n; ++i) {
        if (!(checkEdgeB[i] == networkB.edgeDistribution[i])) {
            edgeB = false;
            break;
        }
    };

    // Overall flag
    bool consistent = nodeA * nodeB * edgeA * edgeB;

    return consistent;
}

// Check convexity of all angles
bool LinkedNetwork::checkConvexity() {
    for (int i = 0; i < networkA.nodes.n; ++i) {
        if (!checkConvexity(i))
            return false;
    }
    return true;
}

// Check convexity by summing angles around node
bool LinkedNetwork::checkConvexity(int id) {
    double angleSum = 0.0;
    // Coordinate vectors, coordination and pbc
    VecF<double> v0(2);
    VecF<double> v1(2);
    VecF<double> v2(2);
    v0[0] = crds[2 * id];
    v0[1] = crds[2 * id + 1];
    int cnd = networkA.nodes[id].netCnxs.n;
    int id1;
    int id2;
    double pbx = networkA.dimensions[0];
    double pby = networkA.dimensions[1];
    double pbrx = networkA.rpb[0];
    double pbry = networkA.rpb[1];
    // Determine vectors to neighbours and sum angles
    for (int i = 0; i < cnd; ++i) {
        int j = (i + 1) % cnd;
        id1 = networkA.nodes[id].netCnxs[i];
        id2 = networkA.nodes[id].netCnxs[j];
        // Periodic vectors to adjacent neighbours
        v1[0] = crds[2 * id1];
        v1[1] = crds[2 * id1 + 1];
        v2[0] = crds[2 * id2];
        v2[1] = crds[2 * id2 + 1];
        v1 -= v0;
        v2 -= v0;
        v1[0] -= pbx * nearbyint(v1[0] * pbrx);
        v1[1] -= pby * nearbyint(v1[1] * pbry);
        v2[0] -= pbx * nearbyint(v2[0] * pbrx);
        v2[1] -= pby * nearbyint(v2[1] * pbry);
        // Angle from dot product
        double n1;
        double n2;
        angleSum += vAngle(v1, v2, n1, n2);
    }
    if (fabs(angleSum - 2 * M_PI) < 1e-12)
        return true;
    else
        return false;
}

// Calculate bond length and angle mean and standard deviation
VecF<double> LinkedNetwork::getOptimisationGeometry(Network network,
                                                    VecF<double> &lenHist,
                                                    VecF<double> &angHist) {

    // Calculate for current configuration
    double x = 0.0;
    double xSq = 0.0;
    double y = 0.0;
    double ySq = 0.0; // len and angle
    int xN = 0;
    int yN = 0; // count
    int cnd;
    VecF<double> v0(2);
    VecF<double> v1(2);
    VecF<double> v2(2);
    double pbx = network.dimensions[0];
    double pby = network.dimensions[1];
    double pbrx = network.rpb[0];
    double pbry = network.rpb[1];
    double lenBinWidth = 4.0 / 10000.0;
    double angBinWidth = 2 * M_PI / 10000.0;
    double bin;
    for (int i = 0; i < network.nodes.n; ++i) {
        cnd = network.nodes[i].netCnxs.n;
        v0[0] = crds[2 * i];
        v0[1] = crds[2 * i + 1];
        for (int j = 0; j < cnd; ++j) {
            int id1 = network.nodes[i].netCnxs[j];
            int id2 = network.nodes[i].netCnxs[(j + 1) % cnd];
            v1[0] = crds[2 * id1];
            v1[1] = crds[2 * id1 + 1];
            v2[0] = crds[2 * id2];
            v2[1] = crds[2 * id2 + 1];
            v1 -= v0;
            v2 -= v0;
            v1[0] -= pbx * nearbyint(v1[0] * pbrx);
            v1[1] -= pby * nearbyint(v1[1] * pbry);
            v2[0] -= pbx * nearbyint(v2[0] * pbrx);
            v2[1] -= pby * nearbyint(v2[1] * pbry);
            double n1;
            double n2;
            double theta = vAngle(v1, v2, n1, n2);
            // Edge lengths avoiding double counting
            if (i < id1) {
                x += n1;
                xSq += n1 * n1;
                xN += 1;
                bin = floor(n1 / lenBinWidth);
                if (bin < lenHist.n)
                    lenHist[bin] += 1.0;
            }
            // Angles
            y += theta;
            ySq += theta * theta;
            yN += 1;
            bin = floor(theta / angBinWidth);
            if (bin < angHist.n)
                angHist[bin] += 1.0;
        }
    }

    // Return current configuration
    VecF<double> optGeom(8);
    optGeom[0] = x;
    optGeom[1] = xSq;
    optGeom[2] = x / xN;
    optGeom[3] = sqrt(xSq / xN - optGeom[2] * optGeom[2]);
    optGeom[4] = y;
    optGeom[5] = ySq;
    optGeom[6] = y / yN;
    optGeom[7] = ySq / yN - optGeom[6] * optGeom[6];
    if (optGeom[7] < 0.0)
        optGeom[7] = 0.0;
    else
        optGeom[7] = sqrt(optGeom[7]);
    return optGeom;
}

// Calculate bond length and angle mean and standard deviation
VecF<double> LinkedNetwork::getOptimisationGeometryTD(VecF<double> &lenHist,
                                                      VecF<double> &angHist) {

    // Calculate for current configuration
    double x = 0.0;
    double xSq = 0.0;
    double y = 0.0;
    double ySq = 0.0; // len and angle
    int xN = 0;
    int yN = 0; // count
    int cnd;
    VecF<double> v0(2);
    VecF<double> v1(2);
    VecF<double> v2(2);
    double pbx = networkT.dimensions[0];
    double pby = networkT.dimensions[1];
    double pbrx = networkT.rpb[0];
    double pbry = networkT.rpb[1];
    double lenBinWidth = 4.0 / 10000.0;
    double angBinWidth = 2 * M_PI / 10000.0;
    double bin;
    for (int i = 0; i < networkT.nodes.n; ++i) {
        cnd = networkT.nodes[i].netCnxs.n;
        v0[0] = crds[2 * i];
        v0[1] = crds[2 * i + 1];
        for (int j = 0; j < cnd; ++j) {
            int id1 = networkT.nodes[i].netCnxs[j];
            int id2 = networkT.nodes[i].netCnxs[(j + 1) % cnd];
            v1[0] = crds[2 * id1];
            v1[1] = crds[2 * id1 + 1];
            v2[0] = crds[2 * id2];
            v2[1] = crds[2 * id2 + 1];
            v1 -= v0;
            v2 -= v0;
            v1[0] -= pbx * nearbyint(v1[0] * pbrx);
            v1[1] -= pby * nearbyint(v1[1] * pbry);
            v2[0] -= pbx * nearbyint(v2[0] * pbrx);
            v2[1] -= pby * nearbyint(v2[1] * pbry);
            double n1;
            double n2;
            double theta = vAngle(v1, v2, n1, n2);
            // Edge lengths avoiding double counting
            if (i < id1) {
                x += n1;
                xSq += n1 * n1;
                xN += 1;
                bin = floor(n1 / lenBinWidth);
                if (bin < lenHist.n)
                    lenHist[bin] += 1.0;
            }
            // Angles
            y += theta;
            ySq += theta * theta;
            yN += 1;
            bin = floor(theta / angBinWidth);
            if (bin < angHist.n)
                angHist[bin] += 1.0;
        }
    }

    // Return current configuration
    VecF<double> optGeom(8);
    optGeom[0] = x;
    optGeom[1] = xSq;
    optGeom[2] = x / xN;
    optGeom[3] = sqrt(xSq / xN - optGeom[2] * optGeom[2]);
    optGeom[4] = y;
    optGeom[5] = ySq;
    optGeom[6] = y / yN;
    optGeom[7] = ySq / yN - optGeom[6] * optGeom[6];
    if (optGeom[7] < 0.0)
        optGeom[7] = 0.0;
    else
        optGeom[7] = sqrt(optGeom[7]);

    return optGeom;
}

// Get sum of areas and squared areas for each ring size
void LinkedNetwork::getRingAreas(VecF<double> &areaSum,
                                 VecF<double> &areaSqSum) {

    // Loop over rings, recentre and apply shoelace formula
    areaSum = 0.0;
    areaSqSum = 0.0;
    double pbx = networkA.dimensions[0];
    double pby = networkA.dimensions[1];
    double pbrx = networkA.rpb[0];
    double pbry = networkA.rpb[1];
    for (int i = 0; i < networkB.nodes.n; ++i) {
        VecR<int> ids = networkB.nodes[i].dualCnxs;
        int ringSize = ids.n;
        VecF<double> xCrds(ringSize);
        VecF<double> yCrds(ringSize);
        for (int j = 0; j < ringSize; ++j) {
            xCrds[j] = crds[2 * ids[j]];
            yCrds[j] = crds[2 * ids[j] + 1];
        }
        xCrds = xCrds - xCrds[0];
        yCrds = yCrds - yCrds[0];
        for (int j = 1; j < ringSize; ++j) {
            xCrds[j] -= pbx * nearbyint(xCrds[j] * pbrx);
            yCrds[j] -= pby * nearbyint(yCrds[j] * pbry);
        }
        double a = 0.0;
        for (int j = 0; j < ringSize - 1; ++j)
            a += xCrds[j] * yCrds[j + 1];
        for (int j = 0; j < ringSize - 1; ++j)
            a -= xCrds[j + 1] * yCrds[j];
        a = 0.5 * fabs(a + xCrds[ringSize - 1] * yCrds[0] -
                       xCrds[0] * yCrds[ringSize - 1]);
        areaSum[ringSize] += a;
        areaSqSum[ringSize] += a * a;
    }
}

// Wrap coordinates of lattice A if periodic
void LinkedNetwork::wrapCoordinates() {
    if (networkT.geometryCode == "2DEtr") {
        HLJ2DP potModel(networkT.dimensions[0], networkT.dimensions[1]);
        potModel.wrap(crds);
    }
    if (networkA.geometryCode == "2DE") {
        HI2DP potModel(networkA.dimensions[0], networkA.dimensions[1]);
        potModel.wrap(crds);
    } else {
        std::ostringstream oss;
        oss << "Unsupported geometry code: " << networkA.geometryCode;
        throw std::runtime_error(oss.str());
    }
}

// Write xyz file format of networks
void LinkedNetwork::writeXYZ(const std::string &prefix) {
    networkA.writeXYZ(prefix + "_A", "O");
    networkB.writeXYZ(prefix + "_B", "N");
    networkT.writeXYZ(prefix + "_T", "Si");
}

// Write networks in format that can be loaded and visualised
void LinkedNetwork::write(const std::string &prefix) {
    networkA.write(prefix + "_A");
    networkB.write(prefix + "_B");
    networkT.write(prefix + "_T");
}
