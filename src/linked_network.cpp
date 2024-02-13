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
    logger->debug("Copying lattice A coordinates");
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
 * @brief Gets the connection type of two connected nodes
 * @param node1Coordination Coordination of the first node
 * @param node2Coordination Coordination of the second node
 * @return The type of the selected connection. 33, 34 or 44.
 * @throw std::runtime_error if the nodes in the random connection have coordinations other than 3 or 4.
 */
int LinkedNetwork::assignValues(int node1Coordination, int node2Coordination) const {
    if (node1Coordination == 3 && node2Coordination == 3) {
        return CNX_TYPE_33;
    } else if (node1Coordination == 4 && node2Coordination == 4) {
        return CNX_TYPE_44;
    } else if ((node1Coordination == 3 && node2Coordination == 4) ||
               (node1Coordination == 4 && node2Coordination == 3)) {
        return CNX_TYPE_43;
    } else {
        std::ostringstream oss;
        oss << "Two nodes have unsupported cooordinations: " << node1Coordination
            << " and " << node2Coordination;
        throw std::runtime_error(oss.str());
    }
}

/**
 * @brief Chooses a random bond in the network and returns the IDs in the bond, the two rings either side and the connection type
 * @param generator Reference to a Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
 * @param selectionType The type of selection to be used. Exponential decay or random.
 * @param logger The logger object.
 * @return A tuple containing the IDs of the two nodes in the bond, the two rings either side and the connection type
 * @throw std::runtime_error if the nodes in the random connection have coordinations other than 3 or 4.
 */
std::tuple<int, int, int, int, int> LinkedNetwork::pickRandomConnection(std::mt19937 &generator, SelectionType selectionType) {
    /*
     * 3-3 coordination connection
     * a,b,c,d,atom4,atom4 are nodes in lattice A
     * ringNode1,ringNode2,ringNode3,ringNode4 are nodes in lattice B
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
    int randNode;
    int randNodeConnection;
    int sharedRingNode1;
    int sharedRingNode2;
    bool pickingAcceptableRing = true;
    while (pickingAcceptableRing) {
        randNode = distribution(generator);
        int randNodeCoordination = networkA.nodes[randNode].netCnxs.n;
        std::uniform_int_distribution<int> randomCnx(0, randNodeCoordination - 1);
        randNodeConnection = networkA.nodes[randNode].netCnxs[randomCnx(generator)];
        int randNodeConnectionCoordination = networkA.nodes[randNodeConnection].netCnxs.n;

        // Use helper function to assign a, b, and cnxType
        cnxType = assignValues(randNodeCoordination, randNodeConnectionCoordination);
        // Two connected nodes should always share two ring nodes.
        if (VecR<int> commonRings = vCommonValues(networkA.nodes[randNode].dualCnxs, networkA.nodes[randNodeConnection].dualCnxs);
            commonRings.n == 2) {
            // Randomly assign ringNode1 and ringNode2 to those two common ring nodes
            std::uniform_int_distribution<int> randomDirection(0, 1);
            int randIndex = randomDirection(generator);
            sharedRingNode1 = commonRings[randIndex];
            sharedRingNode2 = commonRings[1 - randIndex];
        } else {
            wrapCoordinates();
            syncCoordinates();
            write("debug");
            throw std::runtime_error("Selected random connection does not share two ring nodes");
        }
        // Check that the base nodes are not a member of a fixed ring.
        if (!fixedNodes.contains(randNode) && !fixedNodes.contains(randNodeConnection)) {
            // If so, break and return the connection type
            pickingAcceptableRing = false;
        }
    }
    return std::make_tuple(randNode, randNodeConnection, sharedRingNode1, sharedRingNode2, cnxType);
}

/**
 * @brief Generate all ids of nodes in lattices A and B needed for switch move,
 * only for 3/4 coordinate nodes
 * @param cnxType the type of connection to be generated
 * @param switchIDsA the ids of the nodes in lattice A
 * @param switchIDsB the ids of the nodes in lattice B
 * @param switchIDsT the ids of the nodes in lattice T
 * @param baseNode1 the id of the first node in lattice A
 * @param baseNode2 the id of the second node in lattice A
 * @param ringNode1 the id of the first node in lattice B
 * @param ringNode2 the id of the second node in lattice B
 * @return true if the switch move is possible, false otherwise
 */
bool LinkedNetwork::generateSwitchIDs(VecF<int> &switchIDsA, VecF<int> &switchIDsB, VecF<int> &switchIDsT,
                                      int baseNode1, int baseNode2, int ringNode1, int ringNode2, LoggerPtr logger) {

    // lots of error checking to remove any potential pathological cases
    if (baseNode1 == baseNode2 || ringNode1 == ringNode2) {
        logger->warn("Switch move not possible as baseNode1 = baseNode2 or ringNode1 = ringNode2: {} {} {} {}",
                     baseNode1, baseNode2, ringNode1, ringNode2);
        return false;
    }
    /*
     *                7-----8                               7-----8
     *               /       \                              |     |
     *              /         \                      11-----3  2  4-----12
     *      11-----3     2     4-----12                      \   /
     *              \         /                               \ /
     *               \       /                                 1
     *          3     1-----2     4         --->         3     |      4
     *               /       \                                 2
     *              /         \                               /  \
     *      13-----5     1     6-----14                      /    \
     *              \         /                      13-----5  1   6-----14
     *               \       /                              |      |
     *                9-----10                              9------10
     *
     *      Bonds to break       Bonds to Make
     *      1-5, 2-4             1-4, 2-5
     *
     *      Angles to break      Angles to Make
     *      1-5-9, 1-5-13        1-4-8, 1-4-12
     *      2-4-8, 2-4-12        2-5-9, 2-5-13
     *      4-2-1, 4-2-6         4-1-2, 4-1-3
     *      5-1-2, 5-1-3         6-2-1, 6-2-5
     */

    int errorFlag = 0;
    VecR<int> common;
    VecR<int> common1;

    int baseNode5 = findCommonConnection(baseNode1, ringNode1, baseNode2, logger);
    int baseNode6 = findCommonConnection(baseNode2, ringNode1, baseNode1, logger);
    int baseNode3 = findCommonConnection(baseNode1, ringNode2, baseNode2, logger);
    int baseNode4 = findCommonConnection(baseNode2, ringNode2, baseNode1, logger);
    int ringNode3 = findCommonRing(baseNode1, baseNode5, ringNode1, logger);
    int ringNode4 = findCommonRing(baseNode2, baseNode6, ringNode1, logger);

    int baseNode11 = findCommonConnection(baseNode3, ringNode3, baseNode1, logger);
    int baseNode7 = findCommonConnection(baseNode3, ringNode2, baseNode1, logger);
    int baseNode8 = findCommonConnection(baseNode4, ringNode2, baseNode2, logger);
    int baseNode12 = findCommonConnection(baseNode4, ringNode4, baseNode2, logger);
    int baseNode14 = findCommonConnection(baseNode6, ringNode4, baseNode2, logger);
    int baseNode10 = findCommonConnection(baseNode6, ringNode1, baseNode2, logger);
    int baseNode9 = findCommonConnection(baseNode5, ringNode1, baseNode1, logger);
    int baseNode13 = findCommonConnection(baseNode5, ringNode3, baseNode1, logger);

    logger->debug("");
    logger->debug("           {}----{}                                {}------{}", baseNode7, baseNode8, baseNode7, baseNode8);
    logger->debug("          /        \\                                |      |");
    logger->debug("         /          \\                       {}-----{}  {}  {}-----{}", baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
    logger->debug(" {}-----{}    {}    {}-----{}                        \\    /", baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
    logger->debug("         \\          /                                 \\  /");
    logger->debug("          \\        /                                   {}", baseNode1);
    logger->debug("    {}    {}-----{}     {}          --->       {}      |      {}", ringNode3, baseNode1, baseNode2, ringNode4, ringNode3, ringNode4);
    logger->debug("          /       \\                                    {}", baseNode2);
    logger->debug("         /         \\                                  /  \\");
    logger->debug(" {}-----{}    {}    {}-----{}                        /    \\", baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
    logger->debug("         \\         /                        {}-----{}  {}  {}-----{}", baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
    logger->debug("          \\       /                                 |      |");
    logger->debug("          {}----{}                                 {}------{}", baseNode9, baseNode10, baseNode9, baseNode10);
    logger->debug("");

    int alpha = 0;
    int beta = 0;
    int gamma = 0;
    int delta = 0;
    int eta = 0;
    for (int i = 0; i < networkT.nodes[baseNode1].netCnxs.n; ++i) {
        for (int j = 0; j < networkT.nodes[baseNode2].netCnxs.n; ++j) {
            if (networkT.nodes[baseNode1].netCnxs[i] == networkT.nodes[baseNode2].netCnxs[j]) {
                alpha = networkT.nodes[baseNode1].netCnxs[i];
            }
        }
        for (int j = 0; j < networkT.nodes[baseNode5].netCnxs.n; ++j) {
            if (networkT.nodes[baseNode1].netCnxs[i] == networkT.nodes[baseNode5].netCnxs[j]) {
                beta = networkT.nodes[baseNode1].netCnxs[i];
            }
        }
        for (int j = 0; j < networkT.nodes[baseNode3].netCnxs.n; ++j) {
            if (networkT.nodes[baseNode1].netCnxs[i] == networkT.nodes[baseNode3].netCnxs[j]) {
                gamma = networkT.nodes[baseNode1].netCnxs[i];
            }
        }
    }
    for (int i = 0; i < networkT.nodes[baseNode2].netCnxs.n; ++i) {
        for (int j = 0; j < networkT.nodes[baseNode6].netCnxs.n; ++j) {
            if (networkT.nodes[baseNode2].netCnxs[i] == networkT.nodes[baseNode6].netCnxs[j]) {
                delta = networkT.nodes[baseNode2].netCnxs[i];
            }
        }
        for (int j = 0; j < networkT.nodes[baseNode4].netCnxs.n; ++j) {
            if (networkT.nodes[baseNode2].netCnxs[i] == networkT.nodes[baseNode4].netCnxs[j]) {
                eta = networkT.nodes[baseNode2].netCnxs[i];
            }
        }
    }

    // Additional error checking
    if (baseNode5 == baseNode6 || baseNode3 == baseNode4)
        errorFlag = 6; // can simply be triangle edge sharing pair (not an error)

    // Prevent rings having only two or fewer neighbours
    VecR<int> vCnxs = vUnique(networkB.nodes[ringNode2].netCnxs);
    if (vCnxs.n <= 3)
        errorFlag = 10;
    VecR<int> uCnxs = vUnique(networkB.nodes[ringNode1].netCnxs);
    uCnxs.delValue(ringNode2);
    if (uCnxs.n <= 2) {
        for (int i = 0; i < uCnxs.n; ++i) {
            if (uCnxs[i] == ringNode4) {
                errorFlag = 10;
                break;
            }
        }
    }
    if (errorFlag != 0)
        return false;
    // check move will not violate dual connectivity limits
    if (networkB.nodes[ringNode1].netCnxs.n == minBCnxs ||
        networkB.nodes[ringNode2].netCnxs.n == minBCnxs ||
        networkB.nodes[ringNode3].netCnxs.n == maxBCnxs ||
        networkB.nodes[ringNode4].netCnxs.n == maxBCnxs)
        return false;

    switchIDsA = VecF<int>(14);
    switchIDsA[0] = baseNode1;
    switchIDsA[1] = baseNode2;
    switchIDsA[2] = baseNode5;
    switchIDsA[3] = baseNode6;
    switchIDsA[4] = baseNode3;
    switchIDsA[5] = baseNode4;
    switchIDsA[6] = baseNode7;
    switchIDsA[7] = baseNode8;
    switchIDsA[8] = baseNode9;
    switchIDsA[9] = baseNode10;
    switchIDsA[10] = baseNode11;
    switchIDsA[11] = baseNode12;
    switchIDsA[12] = baseNode13;
    switchIDsA[13] = baseNode14;

    switchIDsB = VecF<int>(4);
    switchIDsB[0] = ringNode1;
    switchIDsB[1] = ringNode2;
    switchIDsB[2] = ringNode3;
    switchIDsB[3] = ringNode4;

    switchIDsT = VecF<int>(11);
    switchIDsT[0] = baseNode1;
    switchIDsT[1] = baseNode2;
    switchIDsT[2] = baseNode5;
    switchIDsT[3] = baseNode6;
    switchIDsT[4] = baseNode3;
    switchIDsT[5] = baseNode4;
    switchIDsT[6] = alpha;
    switchIDsT[7] = beta;
    switchIDsT[8] = gamma;
    switchIDsT[9] = delta;
    switchIDsT[10] = eta;
    return true;
}

const int MIN_COORDINATION_NUMBER = 2;
const int NUM_MIX_IDS_A = 6;
const int NUM_MIX_IDS_B = 7;

/**
 * @brief Find a common base node connection between a base node and ring node excluding a given node
 * @param baseNode ID of the base node
 * @param ringNode ID of the ring node
 * @param excludeNode ID of the node to be excluded
 * @return ID of the common connection
 * @throw std::runtime_error if the associated node cannot be found
 */
int LinkedNetwork::findCommonConnection(int baseNode, int ringNode, int excludeNode, LoggerPtr logger) {

    // Find node that shares baseNode and ringNode but is not excludeNode
    int commonConnection = -1;
    VecR<int> commonConnections = vCommonValues(networkA.nodes[baseNode].netCnxs, networkB.nodes[ringNode].dualCnxs);
    commonConnections.delValue(excludeNode);
    if (commonConnections.n == 1)
        commonConnection = commonConnections[0];
    else { // rare high temperature occurrence as a result of 2-cnd nodes giving
           // ring inside ring
        VecR<int> common1(0, commonConnections.n);
        int nCnxs = networkA.nodes[baseNode].netCnxs.n;
        for (int i = 0; i < nCnxs; ++i) {
            if (networkA.nodes[baseNode].netCnxs[i] == excludeNode) {
                int l = networkA.nodes[baseNode].netCnxs[(i + 1) % nCnxs];
                int r = networkA.nodes[baseNode].netCnxs[(i - 1 + nCnxs) % nCnxs];
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
            int nCnxs = networkB.nodes[ringNode].dualCnxs.n;
            for (int i = 0; i < nCnxs; ++i) {
                if (networkB.nodes[ringNode].dualCnxs[i] == excludeNode) {
                    int j;
                    if (networkB.nodes[ringNode].dualCnxs[(i + 1) % nCnxs] == baseNode)
                        j = (i + 2) % nCnxs;
                    else if (networkB.nodes[ringNode].dualCnxs[(i + nCnxs - 1) % nCnxs] == baseNode)
                        j = (i + nCnxs - 2) % nCnxs;
                    if (vContains(common1, networkB.nodes[ringNode].dualCnxs[j])) {
                        commonConnection = networkB.nodes[ringNode].dualCnxs[j];
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
        oss << "Could not find common base node for base node " << baseNode << " and ring node " << ringNode
            << " excluding base node " << excludeNode;
        throw std::runtime_error(oss.str());
    }
    return commonConnection;
}
/**
 * @brief Find a common ring connection between two base nodes that exlcudes a given node
 * @param baseNode1 ID of the first base node
 * @param baseNode2 ID of the second base node
 * @param excludeNode ID of the ring node to be excluded
 * @return ID of the common ring connection
 * @throw std::runtime_error if the associated node cannot be found
 */
int LinkedNetwork::findCommonRing(int baseNode1, int baseNode2, int excludeNode, LoggerPtr logger) {

    // Find node that shares baseNode1 and baseNode2 but is not excludeNode
    int associated = -1;
    VecR<int> common;
    common = vCommonValues(networkA.nodes[baseNode1].dualCnxs, networkA.nodes[baseNode2].dualCnxs);
    common.delValue(excludeNode);
    if (common.n == 1)
        associated = common[0];
    else { // rare case with periodic interactions
        VecR<int> common1(0, common.n);
        for (int i = 0; i < common.n; ++i) {
            if (vContains(networkB.nodes[common[i]].netCnxs, excludeNode))
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
                    if ((ring[j] == baseNode1 && ring[k] == baseNode2) || (ring[j] == baseNode2 && ring[k] == baseNode1)) {
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
        std::ostringstream oss;
        oss << "Could not find common ring node for base node " << baseNode1 << " and base node " << baseNode2
            << " excluding ring node " << excludeNode;
        throw std::runtime_error(oss.str());
    }
    return associated;
}

// Switch connectivities in lattice between 2x3 coordinate nodes
void LinkedNetwork::switchCnx33(VecF<int> switchIDsA, VecF<int> switchIDsB,
                                VecF<int> switchIDsT, LoggerPtr logger) {

    /*
     *                7-----8                               7-----8
     *               /       \                              |     |
     *              /         \                      11-----3  2  4-----12
     *      11-----3     2     4-----12                      \   /
     *              \         /                               \ /
     *               \       /                                 1
     *          3     1-----2     4         --->         3     |      4
     *               /       \                                 2
     *              /         \                               /  \
     *      13-----5     1     6-----14                      /    \
     *              \         /                      13-----5  1   6-----14
     *               \       /                              |      |
     *                9-----10                              9------10
     *
     *      Bonds to break       Bonds to Make
     *      1-5, 2-4             1-4, 2-5
     *
     *      Angles to break      Angles to Make
     *      1-5-9, 1-5-13        1-4-8, 1-4-12
     *      2-4-8, 2-4-12        2-5-9, 2-5-13
     *      4-2-1, 4-2-6         4-1-2, 4-1-3
     *      5-1-2, 5-1-3         6-2-1, 6-2-5
     */

    /**
     * 0 -> a -> baseNode1
     * 1 -> b -> baseNode2
     * 2 -> c -> baseNode5
     * 3 -> d -> baseNode6
     * 4 -> e -> baseNode3
     * 5 -> f -> baseNode4
     *
     * 1 -> u -> ringNode1
     * 2 -> v -> ringNode2
     * 3 -> w -> ringNode3
     * 4 -> x -> ringNode4
     */

    int atom1 = switchIDsA[0];
    int atom2 = switchIDsA[1];
    int atom5 = switchIDsA[2];
    int atom4 = switchIDsA[5];

    int ringNode1 = switchIDsB[0];
    int ringNode2 = switchIDsB[1];
    int ringNode3 = switchIDsB[2];
    int ringNode4 = switchIDsB[3];

    int beta = switchIDsT[7];
    int gamma = switchIDsT[8];
    int delta = switchIDsT[9];
    int eta = switchIDsT[10];
    // Apply changes to descriptors due to breaking connections
    // For network A node distribution and edge distribution will remain unchanged

    // For network B node and edge distribution will change
    int nu = networkB.nodes[ringNode1].netCnxs.n;
    int nv = networkB.nodes[ringNode2].netCnxs.n;
    int nw = networkB.nodes[ringNode3].netCnxs.n;
    int nx = networkB.nodes[ringNode4].netCnxs.n;
    --networkB.nodeDistribution[nu];
    --networkB.nodeDistribution[nv];
    --networkB.nodeDistribution[nw];
    --networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[ringNode1].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nu][nCnx];
        if (id != ringNode2 && id != ringNode3 && id != ringNode4)
            --networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[ringNode2].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nv][nCnx];
        if (id != ringNode1 && id != ringNode3 && id != ringNode4)
            --networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[ringNode3].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nw][nCnx];
        if (id != ringNode1 && id != ringNode2 && id != ringNode4)
            --networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[ringNode4].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        --networkB.edgeDistribution[nx][nCnx];
        if (id != ringNode1 && id != ringNode2 && id != ringNode3)
            --networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }

    logger->debug("Switching A-A connections");
    // A-A connectivities

    networkA.nodes[atom1].netCnxs.replaceValue(atom5, atom4);
    networkA.nodes[atom2].netCnxs.replaceValue(atom4, atom5);
    networkA.nodes[atom4].netCnxs.replaceValue(atom2, atom1);
    networkA.nodes[atom5].netCnxs.replaceValue(atom1, atom2);

    logger->debug("Switching A-B connections");
    // A-B connectvities
    networkA.nodes[atom1].dualCnxs.replaceValue(ringNode1, ringNode4);
    networkA.nodes[atom2].dualCnxs.replaceValue(ringNode2, ringNode3);

    logger->debug("Switching B-B connections");
    // B-B connectivities
    networkB.nodes[ringNode1].netCnxs.delValue(ringNode2);
    networkB.nodes[ringNode2].netCnxs.delValue(ringNode1);
    networkB.nodes[ringNode3].netCnxs.insertValue(ringNode4, ringNode1, ringNode2);
    networkB.nodes[ringNode4].netCnxs.insertValue(ringNode3, ringNode1, ringNode2);

    logger->debug("Switching B-A connections");
    // B-A connectivities
    networkB.nodes[ringNode1].dualCnxs.delValue(atom1);
    networkB.nodes[ringNode2].dualCnxs.delValue(atom2);
    networkB.nodes[ringNode3].dualCnxs.insertValue(atom2, atom1, atom5);
    networkB.nodes[ringNode4].dualCnxs.insertValue(atom1, atom2, atom4);

    // Apply changes to descriptors due to making connections
    // Network B
    nu = networkB.nodes[ringNode1].netCnxs.n;
    nv = networkB.nodes[ringNode2].netCnxs.n;
    nw = networkB.nodes[ringNode3].netCnxs.n;
    nx = networkB.nodes[ringNode4].netCnxs.n;
    ++networkB.nodeDistribution[nu];
    ++networkB.nodeDistribution[nv];
    ++networkB.nodeDistribution[nw];
    ++networkB.nodeDistribution[nx];
    for (int i = 0; i < nu; ++i) {
        int id = networkB.nodes[ringNode1].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nu][nCnx];
        if (id != ringNode2 && id != ringNode3 && id != ringNode4)
            ++networkB.edgeDistribution[nCnx][nu]; // prevent double counting
    }
    for (int i = 0; i < nv; ++i) {
        int id = networkB.nodes[ringNode2].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nv][nCnx];
        if (id != ringNode1 && id != ringNode3 && id != ringNode4)
            ++networkB.edgeDistribution[nCnx][nv]; // prevent double counting
    }
    for (int i = 0; i < nw; ++i) {
        int id = networkB.nodes[ringNode3].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nw][nCnx];
        if (id != ringNode1 && id != ringNode2 && id != ringNode4)
            ++networkB.edgeDistribution[nCnx][nw]; // prevent double counting
    }
    for (int i = 0; i < nx; ++i) {
        int id = networkB.nodes[ringNode4].netCnxs[i];
        int nCnx = networkB.nodes[id].netCnxs.n;
        ++networkB.edgeDistribution[nx][nCnx];
        if (id != ringNode1 && id != ringNode2 && id != ringNode3)
            ++networkB.edgeDistribution[nCnx][nx]; // prevent double counting
    }

    logger->debug("Switching network T...");
    networkT.nodes[atom1].netCnxs.replaceValue(gamma, delta);
    networkT.nodes[atom2].netCnxs.replaceValue(delta, gamma);
    networkT.nodes[delta].netCnxs.replaceValue(atom2, atom1);
    networkT.nodes[gamma].netCnxs.replaceValue(atom1, atom2);
    networkT.nodes[gamma].netCnxs.replaceValue(beta, eta);
    networkT.nodes[beta].netCnxs.replaceValue(gamma, delta);
    networkT.nodes[delta].netCnxs.replaceValue(eta, beta);
    networkT.nodes[eta].netCnxs.replaceValue(delta, gamma);
    networkT.nodes[atom1].dualCnxs.replaceValue(ringNode2, ringNode4);
    networkT.nodes[atom2].dualCnxs.replaceValue(ringNode1, ringNode3);
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
bool LinkedNetwork::convexRearrangement(VecF<int> switchIdsA, LoggerPtr logger) {

    bool convex;
    // Unpack nodes
    int baseNode1 = switchIdsA[0];
    int baseNode2 = switchIdsA[1];
    int baseNode5 = switchIdsA[2];
    int baseNode6 = switchIdsA[3];
    int baseNode3 = switchIdsA[4];
    int baseNode4 = switchIdsA[5];

    // Maintains which nodes are convex
    VecF<bool> convexNodes;
    convexNodes = VecF<bool>(6);

    // Initial guess places baseNode1 at the centre of cd, baseNode2 at the centre of ef
    VecF<double> va(2);
    VecF<double> vb(2);
    VecF<double> vc(2);
    VecF<double> vd(2);
    VecF<double> ve(2);
    VecF<double> vf(2);
    va[0] = crds[2 * baseNode1];
    va[1] = crds[2 * baseNode1 + 1];
    vb[0] = crds[2 * baseNode2];
    vb[1] = crds[2 * baseNode2 + 1];
    vc[0] = crds[2 * baseNode5];
    vc[1] = crds[2 * baseNode5 + 1];
    vd[0] = crds[2 * baseNode6];
    vd[1] = crds[2 * baseNode6 + 1];
    ve[0] = crds[2 * baseNode3];
    ve[1] = crds[2 * baseNode3 + 1];
    vf[0] = crds[2 * baseNode4];
    vf[1] = crds[2 * baseNode4 + 1];
    VecF<double> vce(2);
    VecF<double> vdf(2);
    VecF<double> vcd(2);
    VecF<double> vef(2);
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
    crds[2 * baseNode1] = va[0];
    crds[2 * baseNode1 + 1] = va[1];
    crds[2 * baseNode2] = vb[0];
    crds[2 * baseNode2 + 1] = vb[1];
    for (int i = 0; i < 6; ++i)
        convexNodes[i] = checkConvexity(switchIdsA[i]);
    convex = (convexNodes == true);

    // Guess move baseNode1,baseNode2 towards each other
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
            crds[2 * baseNode1] = va[0];
            crds[2 * baseNode1 + 1] = va[1];
            crds[2 * baseNode2] = vb[0];
            crds[2 * baseNode2 + 1] = vb[1];
            for (int k = 0; k < 6; ++k)
                convexNodes[k] = checkConvexity(switchIdsA[k]);
            convex = (convexNodes == true);
            if (convex)
                break;
        }
    }
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
    logger->debug("Finding move...");

    VecF<int> switchIDsA;
    VecF<int> switchIDsB;
    VecF<int> switchIDsT;
    bool foundValidMove = false;
    int baseNode1;
    int baseNode2;
    int ringNode1;
    int ringNode2;
    int cnxType;
    for (int i = 0; i < networkA.nodes.n * networkA.nodes.n; ++i) {
        std::tuple<int, int, int, int, int> result;
        if (mcWeighting == "weighted")
            result = pickRandomConnection(mtGen, SelectionType::EXPONENTIAL_DECAY);
        else
            result = pickRandomConnection(mtGen, SelectionType::RANDOM);

        std::tie(baseNode1, baseNode2, ringNode1, ringNode2, cnxType) = result;
        logger->debug("Picked base nodes: {} {} and ring nodes: {} {}", baseNode1, baseNode2, ringNode1, ringNode2);
        foundValidMove = generateSwitchIDs(switchIDsA, switchIDsB, switchIDsT, baseNode1, baseNode2, ringNode1, ringNode2, logger);
        if (foundValidMove)
            break;
    }

    // Now baseNode1, baseNode2, ringNode1, ringNode2, and cnxType are available here
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

    VecF<Node> saveNodesA(switchIDsA.n);
    VecF<Node> saveNodesB(switchIDsB.n);
    VecF<Node> saveNodesT(switchIDsT.n);
    for (int i = 0; i < saveNodesA.n; ++i)
        saveNodesA[i] = networkA.nodes[switchIDsA[i]];

    for (int i = 0; i < saveNodesB.n; ++i)
        saveNodesB[i] = networkB.nodes[switchIDsB[i]];
    for (int i = 0; i < saveNodesT.n; ++i)
        saveNodesT[i] = networkT.nodes[switchIDsT[i]];

    // Switch and geometry optimise
    VecF<int> optStatus_networkA(2);
    VecF<int> optStatus_SimpleGraphene(2);
    VecF<int> optStatus_TersoffGraphene(2);
    VecF<int> optStatus_TriangleRaft(2);
    VecF<int> optStatus_Bilayer(2);
    VecF<int> optStauts_BN(2);

    // works for network version of lammps objects
    logger->debug("Switching NetMC Network...");
    switchCnx33(switchIDsA, switchIDsB, switchIDsT, logger);
    // works for lammps objects

    logger->debug("Switching LAMMPS Networks...");
    if (isSimpleGrapheneEnabled) {
        logger->debug("Switching Simple Graphene");
        SimpleGraphene.switchGraphene(switchIDsA, logger);
    }
    if (isTriangleRaftEnabled) {
        logger->debug("Switching Triangle Raft");
        Triangle_Raft.switchTriangleRaft(switchIDsA, switchIDsT, logger);
    }
    if (isBilayerEnabled)
        Bilayer.switchBilayer(switchIDsA, switchIDsT, logger);
    if (isBNEnabled) {
        logger->debug("Switching BN");
        BN.switchGraphene(switchIDsA, logger);
    }
    // Rearrange nodes after switch
    bool geometryOK = true;
    geometryOK = checkThreeRingEdges(ringNode1);
    if (geometryOK)
        geometryOK = checkThreeRingEdges(ringNode2);
    logger->debug("1");
    if (geometryOK) {
        logger->debug("2");
        if (isMaintainConvexityEnabled) {
            logger->debug("3");
            geometryOK = convexRearrangement(switchIDsA, logger);
            for (int i = 0; i < 6; ++i) {
                geometryOK = checkConvexity(switchIDsA[i]);
                if (!geometryOK)
                    break;
            }
        }
    } else {
        optStatus_SimpleGraphene = VecF<int>(3);
    }
    logger->debug("4");
    if (!geometryOK)
        optStatus_SimpleGraphene[0] = 4;
    int nSi;

    // Geometry optimisation of local region
    logger->debug("Optimising geometry...");
    if (geometryOK) {
        logger->debug("Geometry OK");
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
        }

        for (int i = 0; i < 6; ++i) {
            if (isSimpleGrapheneEnabled && isTersoffGrapheneEnabled) {
                localCrdsTersoff[2 * switchIDsA[i]] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i]] * CScaling;
                localCrdsTersoff[2 * switchIDsA[i] + 1] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i] + 1] * CScaling;
            }
            if (isSimpleGrapheneEnabled && isBilayerEnabled) {
                localCrdsBilayer[3 * (switchIDsA[i] * 2)] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIDsA[i] * 2) + 1] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i] + 1] * SiScaling;
                localCrdsBilayer[3 * (switchIDsA[i] * 2 + 1)] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIDsA[i] * 2 + 1) + 1] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i] + 1] * SiScaling;
                localCrdsBilayer[3 * (switchIDsA[i] + nSi)] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i]] * SiScaling;
                localCrdsBilayer[3 * (switchIDsA[i] + nSi) + 1] =
                    localCrdsSimpleGraphene[2 * switchIDsA[i] + 1] * SiScaling;
            }
        }
        if (isTersoffGrapheneEnabled)
            TersoffGraphene.pushCrds(2, localCrdsTersoff);
        if (isBilayerEnabled)
            Bilayer.pushCrds(3, localCrdsBilayer);

        if (isOpenMPIEnabled) {
            optStatus_networkA = localGeometryOptimisation(
                baseNode1, baseNode2, goptParamsA[1], false, isMaintainConvexityEnabled, logger);
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
                baseNode1, baseNode2, goptParamsA[1], 0, isMaintainConvexityEnabled, logger);
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
            logger->info("Rejected Move, Ei = {} Ef = {}", saveEnergySimpleGraphene, SimpleGrapheneEnergy);
        else if (mcRoutine == 2)
            logger->info("Rejected Move, Ei = {} Ef = {}", saveEnergyTriangleRaft, TriangleRaftEnergy);
        else if (mcRoutine == 5)
            logger->info("Rejected Move, Ei = {} Ef = {}", saveEnergyBN, BNEnergy);

        crds = saveCrds;
        networkA.nodeDistribution = saveNodeDistA;
        networkA.edgeDistribution = saveEdgeDistA;
        networkB.nodeDistribution = saveNodeDistB;
        networkB.edgeDistribution = saveEdgeDistB;

        for (int i = 0; i < saveNodesA.n; ++i)
            networkA.nodes[switchIDsA[i]] = saveNodesA[i];
        for (int i = 0; i < saveNodesB.n; ++i)
            networkB.nodes[switchIDsB[i]] = saveNodesB[i];
        for (int i = 0; i < saveNodesT.n; ++i)
            networkT.nodes[switchIDsT[i]] = saveNodesT[i];

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
            SimpleGraphene.revertGraphene(switchIDsA, logger);
            SimpleGraphene.globalPotentialMinimisation();
        }
        if (isTriangleRaftEnabled) {
            Triangle_Raft.revertTriangleRaft(switchIDsA, switchIDsT, logger);
            Triangle_Raft.globalPotentialMinimisation();
        }
        if (isBilayerEnabled) {
            Bilayer.revertBilayer(switchIDsA, switchIDsT, logger);
            Bilayer.globalPotentialMinimisation();
        }
        if (isBNEnabled) {
            BN.revertGraphene(switchIDsA, logger);
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
    int atom1;
    int atom2;
    int ringNode1;
    int ringNode2;
    VecF<int> switchIDsA;
    VecF<int> switchIDsB;
    VecF<int> switchIDsT;
    bool foundValidMove = false;
    int cnxType;
    for (int i = 0; i < networkA.nodes.n * networkA.nodes.n;
         ++i) { // catch in case cannot find any valid moves
        std::tuple<int, int, int, int, int> result;
        result = pickRandomConnection(mtGen, SelectionType::RANDOM);
        std::tie(atom1, atom2, ringNode1, ringNode2, cnxType) = result;
        foundValidMove = generateSwitchIDs(switchIDsA, switchIDsB,
                                           switchIDsT, atom1, atom2, ringNode1, ringNode2, logger);
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
    VecF<Node> saveNodesA(switchIDsA.n);
    VecF<Node> saveNodesB(switchIDsB.n);
    VecF<Node> saveNodesT(switchIDsT.n);
    for (int i = 0; i < saveNodesA.n; ++i)
        saveNodesA[i] = networkA.nodes[switchIDsA[i]];
    for (int i = 0; i < saveNodesB.n; ++i)
        saveNodesB[i] = networkB.nodes[switchIDsB[i]];
    for (int i = 0; i < saveNodesT.n; ++i)
        saveNodesT[i] = networkT.nodes[switchIDsT[i]];

    // Switch and geometry optimise
    logger->info("Switching...");
    VecF<int> optStatus(3);
    switchCnx33(switchIDsA, switchIDsB, switchIDsT, logger);

    // Rearrange nodes after switch
    bool geometryOK = true;
    geometryOK = checkThreeRingEdges(ringNode1);
    if (geometryOK)
        geometryOK = checkThreeRingEdges(ringNode2);
    if (geometryOK) {
        if (isMaintainConvexityEnabled) {
            geometryOK = convexRearrangement(switchIDsA, logger);
            for (int i = 0; i < 6; ++i) {
                geometryOK = checkConvexity(switchIDsA[i]);
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
            networkA.nodes[switchIDsA[i]] = saveNodesA[i];
        for (int i = 0; i < saveNodesB.n; ++i)
            networkB.nodes[switchIDsB[i]] = saveNodesB[i];
        for (int i = 0; i < saveNodesT.n; ++i)
            networkT.nodes[switchIDsT[i]] = saveNodesT[i];
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
        VecF<double> ringNode4(networkB.nodes[i].dualCnxs.n);
        VecF<double> y(networkB.nodes[i].dualCnxs.n);
        for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
            ringNode4[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
            y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
        }
        VecF<double> origin(2);
        origin[0] = ringNode4[0];
        origin[1] = y[0];
        ringNode4 -= origin[0];
        y -= origin[1];
        for (int j = 0; j < ringNode4.n; ++j)
            ringNode4[j] -= networkB.dimensions[0] * nearbyint(ringNode4[j] * networkB.rpb[0]);
        for (int j = 0; j < y.n; ++j)
            y[j] -= networkB.dimensions[1] * nearbyint(y[j] * networkB.rpb[1]);
        VecF<double> c(2);
        c[0] = origin[0] + vMean(ringNode4);
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
        VecF<double> ringNode4(networkB.nodes[i].dualCnxs.n);
        VecF<double> y(networkB.nodes[i].dualCnxs.n);
        for (int j = 0; j < networkB.nodes[i].dualCnxs.n; ++j) {
            ringNode4[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[0];
            y[j] = networkA.nodes[networkB.nodes[i].dualCnxs[j]].crd[1];
        }
        VecF<double> origin(2);
        origin[0] = ringNode4[0];
        origin[1] = y[0];
        ringNode4 -= origin[0];
        y -= origin[1];
        for (int j = 0; j < ringNode4.n; ++j)
            ringNode4[j] -= networkB.dimensions[0] * nearbyint(ringNode4[j] * networkB.rpb[0]);
        for (int j = 0; j < y.n; ++j)
            y[j] -= networkB.dimensions[1] * nearbyint(y[j] * networkB.rpb[1]);
        VecF<double> c(2);
        c[0] = origin[0] + vMean(ringNode4);
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

// Check if the polygon formed by the connections of a node is convex
/**
bool LinkedNetwork::checkConvexity(int nodeId) {
    double totalAngle = 0.0;
    VecF<double> nodeCoords(2);
    VecF<double> firstNeighborCoords(2);
    VecF<double> secondNeighborCoords(2);

    // Get the coordinates of the node
    nodeCoords[0] = crds[2 * nodeId];
    nodeCoords[1] = crds[2 * nodeId + 1];

    // Get the number of connections of the node
    int numConnections = networkA.nodes[nodeId].netCnxs.n;

    // Get the dimensions of the network and the reciprocal of the dimensions
    double networkWidth = networkA.dimensions[0];
    double networkHeight = networkA.dimensions[1];
    double reciprocalWidth = networkA.rpb[0];
    double reciprocalHeight = networkA.rpb[1];

    // Iterate over each connection of the node
    for (int i = 0; i < numConnections; ++i) {
        int nextIndex = (i + 1) % numConnections;
        int firstNeighborId = networkA.nodes[nodeId].netCnxs[i];
        int secondNeighborId = networkA.nodes[nodeId].netCnxs[nextIndex];

        // Get the coordinates of the first and second neighbors
        firstNeighborCoords[0] = crds[2 * firstNeighborId];
        firstNeighborCoords[1] = crds[2 * firstNeighborId + 1];
        secondNeighborCoords[0] = crds[2 * secondNeighborId];
        secondNeighborCoords[1] = crds[2 * secondNeighborId + 1];

        // Adjust the coordinates for periodic boundary conditions
        firstNeighborCoords -= nodeCoords;
        secondNeighborCoords -= nodeCoords;
        firstNeighborCoords[0] -= networkWidth * nearbyint(firstNeighborCoords[0] * reciprocalWidth);
        firstNeighborCoords[1] -= networkHeight * nearbyint(firstNeighborCoords[1] * reciprocalHeight);
        secondNeighborCoords[0] -= networkWidth * nearbyint(secondNeighborCoords[0] * reciprocalWidth);
        secondNeighborCoords[1] -= networkHeight * nearbyint(secondNeighborCoords[1] * reciprocalHeight);

        // Calculate the angle between the first and second neighbors and add it to the total angle
        double n1;
        double n2;
        totalAngle += vAngle(firstNeighborCoords, secondNeighborCoords, n1, n2);
    }

    // If the total angle is 2 * PI (360 degrees), the polygon is convex
    if (fabs(totalAngle - 2 * M_PI) < 1e-12)
        return true;
    else
        return false;
}
*/

// Calculate bond length and angle mean and standard deviation
VecF<double> LinkedNetwork::getOptimisationGeometry(Network network,
                                                    VecF<double> &lenHist,
                                                    VecF<double> &angHist) {

    // Calculate for current configuration
    double ringNode4 = 0.0;
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
                ringNode4 += n1;
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
    optGeom[0] = ringNode4;
    optGeom[1] = xSq;
    optGeom[2] = ringNode4 / xN;
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
    double ringNode4 = 0.0;
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
                ringNode4 += n1;
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
    optGeom[0] = ringNode4;
    optGeom[1] = xSq;
    optGeom[2] = ringNode4 / xN;
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
        double area = 0.0;
        for (int j = 0; j < ringSize - 1; ++j)
            area += xCrds[j] * yCrds[j + 1];
        for (int j = 0; j < ringSize - 1; ++j)
            area -= xCrds[j + 1] * yCrds[j];
        area = 0.5 * fabs(area + xCrds[ringSize - 1] * yCrds[0] -
                          xCrds[0] * yCrds[ringSize - 1]);
        areaSum[ringSize] += area;
        areaSqSum[ringSize] += area * area;
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
