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
    currentCoords.resize(2 * networkA.nodes.n);
    networkA.getCoords(currentCoords);
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

    if (inputData.isSimpleGrapheneEnabled)
        lammpsNetwork = LammpsObject("Si", inputData.inputFolder, logger);
    else if (inputData.isTersoffGrapheneEnabled)
        lammpsNetwork = LammpsObject("C", inputData.inputFolder, logger);
    else if (inputData.isTriangleRaftEnabled)
        lammpsNetwork = LammpsObject("Si2O3", inputData.inputFolder, logger);
    else if (inputData.isBilayerEnabled)
        lammpsNetwork = LammpsObject("SiO2", inputData.inputFolder, logger);
    else if (inputData.isBNEnabled)
        lammpsNetwork = LammpsObject("BN", inputData.inputFolder, logger);
    else
        throw std::runtime_error("You must enable at least one network type");
    if (inputData.isFixRingsEnabled) {
        findFixedRings(inputData.isFixRingsEnabled, prefix, logger);
        findFixedNodes();
    }

    isMaintainConvexityEnabled = inputData.isMaintainConvexityEnabled;
    energy = lammpsNetwork.getPotentialEnergy();
    currentCoords.resize(2 * networkA.nodes.n);

    networkA.getCoords(currentCoords);

    mc = Metropolis(inputData.randomSeed, pow(10, inputData.startTemperature), energy);
    mtGen.seed(inputData.randomSeed);

    isOpenMPIEnabled = inputData.isOpenMPIEnabled;
    mcWeighting = inputData.randomOrWeighted;
    isSimpleGrapheneEnabled = inputData.isSimpleGrapheneEnabled;
    isTersoffGrapheneEnabled = inputData.isTersoffGrapheneEnabled;
    isTriangleRaftEnabled = inputData.isTriangleRaftEnabled;
    isBilayerEnabled = inputData.isBilayerEnabled;
    isBNEnabled = inputData.isBNEnabled;
    mcRoutine = inputData.selectedMinimisationProtocol;
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

// Single monte carlo switching move
void LinkedNetwork::monteCarloSwitchMoveLAMMPS(LoggerPtr logger) {

    /* Single MC switch move
     * 1) select random connection
     * 2) switch connection
     * 3) optimise and evaluate energy
     * 4) accept or reject */
    logger->debug("Finding move...");

    VecF<int> switchIDsA;
    VecF<int> switchIDsB;
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
        foundValidMove = generateSwitchIDs(switchIDsA, switchIDsB, baseNode1, baseNode2, ringNode1, ringNode2, logger);
        if (foundValidMove)
            break;
    }

    // Now baseNode1, baseNode2, ringNode1, ringNode2, and cnxType are available here
    if (!foundValidMove) {
        logger->error("Cannot find any valid switch moves");
        throw std::runtime_error("Cannot find any valid switch moves");
    }
    numSwitches++;
    // Save current state
    double initialEnergy = energy;

    logger->debug("Saving initial coordinates...");
    std::vector<double> initialCoords = currentCoords;

    VecF<int> saveNodeDistA = networkA.nodeDistribution;
    VecF<int> saveNodeDistB = networkB.nodeDistribution;

    VecF<VecF<int>> saveEdgeDistA = networkA.edgeDistribution;
    VecF<VecF<int>> saveEdgeDistB = networkB.edgeDistribution;

    VecF<Node> saveNodesA(switchIDsA.n);
    VecF<Node> saveNodesB(switchIDsB.n);
    for (int i = 0; i < saveNodesA.n; ++i)
        saveNodesA[i] = networkA.nodes[switchIDsA[i]];

    for (int i = 0; i < saveNodesB.n; ++i)
        saveNodesB[i] = networkB.nodes[switchIDsB[i]];

    // Switch and geometry optimise

    logger->debug("Switching NetMC Network...");
    switchCnx33(switchIDsA, switchIDsB, logger);

    logger->debug("Switching LAMMPS Networks...");
    if (isSimpleGrapheneEnabled) {
        logger->debug("Switching Simple Graphene");
        lammpsNetwork.switchGraphene(switchIDsA, logger);
    }
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
    }
    double finalEnergy;
    // Geometry optimisation of local region
    logger->debug("Optimising geometry...");
    if (geometryOK) {
        logger->debug("Geometry OK");
        lammpsNetwork.minimiseNetwork();
        logger->debug("A");
        finalEnergy = lammpsNetwork.getPotentialEnergy();
        logger->debug("B");
    } else {
        logger->debug("Geometry not OK, setting final energy to infinity.");
        finalEnergy = std::numeric_limits<double>::infinity();
    }
    logger->info("Accepting or rejecting...");
    bool isAccepted = mc.acceptanceCriterion(finalEnergy, initialEnergy, 1.0);
    if (isAccepted) {
        logger->info("Accepted Move, Ei = {} Ef = {}", initialEnergy, finalEnergy);
        logger->info("Syncing coordinates...");
        currentCoords = lammpsNetwork.getCoords(2);
        wrapCoords(currentCoords);
        pushCoords(currentCoords);
        energy = finalEnergy;
        numAcceptedSwitches++;
    } else {
        logger->info("Rejected Move, Ei = {} Ef = {}", initialEnergy, finalEnergy);
        lammpsNetwork.revertGraphene(switchIDsA, logger);
        lammpsNetwork.setCoords(initialCoords, 2);

        networkA.nodeDistribution = saveNodeDistA;
        networkA.edgeDistribution = saveEdgeDistA;
        networkB.nodeDistribution = saveNodeDistB;
        networkB.edgeDistribution = saveEdgeDistB;

        for (int i = 0; i < saveNodesA.n; ++i)
            networkA.nodes[switchIDsA[i]] = saveNodesA[i];
        for (int i = 0; i < saveNodesB.n; ++i)
            networkB.nodes[switchIDsB[i]] = saveNodesB[i];
    }
}

void LinkedNetwork::showCoords(std::vector<double> &coords, LoggerPtr logger) {
    for (int i = 0; i < coords.size(); i += 2) {
        logger->debug("{}) {} {}", i, coords[i], coords[i + 1]);
    }
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
 * @param baseNode1 the id of the first node in lattice A
 * @param baseNode2 the id of the second node in lattice A
 * @param ringNode1 the id of the first node in lattice B
 * @param ringNode2 the id of the second node in lattice B
 * @return true if the switch move is possible, false otherwise
 */
bool LinkedNetwork::generateSwitchIDs(VecF<int> &switchIDsA, VecF<int> &switchIDsB,
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
        std::ostringstream oss;
        oss << "Could not find common ring node for base node " << baseNode1 << " and base node " << baseNode2
            << " excluding ring node " << excludeNode;
        throw std::runtime_error(oss.str());
    }
    return associated;
}

// Switch connectivities in lattice between 2x3 coordinate nodes
void LinkedNetwork::switchCnx33(VecF<int> switchIDsA, VecF<int> switchIDsB, LoggerPtr logger) {

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

    // logger->debug("Switching A-A connections");
    // A-A connectivities

    networkA.nodes[atom1].netCnxs.replaceValue(atom5, atom4);
    networkA.nodes[atom2].netCnxs.replaceValue(atom4, atom5);
    networkA.nodes[atom4].netCnxs.replaceValue(atom2, atom1);
    networkA.nodes[atom5].netCnxs.replaceValue(atom1, atom2);

    // logger->debug("Switching A-B connections");
    // A-B connectvities
    networkA.nodes[atom1].dualCnxs.replaceValue(ringNode1, ringNode4);
    networkA.nodes[atom2].dualCnxs.replaceValue(ringNode2, ringNode3);

    // logger->debug("Switching B-B connections");
    // B-B connectivities
    networkB.nodes[ringNode1].netCnxs.delValue(ringNode2);
    networkB.nodes[ringNode2].netCnxs.delValue(ringNode1);
    networkB.nodes[ringNode3].netCnxs.insertValue(ringNode4, ringNode1, ringNode2);
    networkB.nodes[ringNode4].netCnxs.insertValue(ringNode3, ringNode1, ringNode2);

    // logger->debug("Switching B-A connections");
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
    va[0] = currentCoords[2 * baseNode1];
    va[1] = currentCoords[2 * baseNode1 + 1];
    vb[0] = currentCoords[2 * baseNode2];
    vb[1] = currentCoords[2 * baseNode2 + 1];
    vc[0] = currentCoords[2 * baseNode5];
    vc[1] = currentCoords[2 * baseNode5 + 1];
    vd[0] = currentCoords[2 * baseNode6];
    vd[1] = currentCoords[2 * baseNode6 + 1];
    ve[0] = currentCoords[2 * baseNode3];
    ve[1] = currentCoords[2 * baseNode3 + 1];
    vf[0] = currentCoords[2 * baseNode4];
    vf[1] = currentCoords[2 * baseNode4 + 1];
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
    currentCoords[2 * baseNode1] = va[0];
    currentCoords[2 * baseNode1 + 1] = va[1];
    currentCoords[2 * baseNode2] = vb[0];
    currentCoords[2 * baseNode2 + 1] = vb[1];
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
            currentCoords[2 * baseNode1] = va[0];
            currentCoords[2 * baseNode1 + 1] = va[1];
            currentCoords[2 * baseNode2] = vb[0];
            currentCoords[2 * baseNode2 + 1] = vb[1];
            for (int k = 0; k < 6; ++k)
                convexNodes[k] = checkConvexity(switchIdsA[k]);
            convex = (convexNodes == true);
            if (convex)
                break;
        }
    }
    return convex;
}
void LinkedNetwork::pushCoords(std::vector<double> coords) {

    // Sync A coordinates
    for (int i = 0; i < networkA.nodes.n; ++i) {
        networkA.nodes[i].crd[0] = coords[2 * i];
        networkA.nodes[i].crd[1] = coords[2 * i + 1];
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
    v0[0] = currentCoords[2 * id];
    v0[1] = currentCoords[2 * id + 1];
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
        v1[0] = currentCoords[2 * id1];
        v1[1] = currentCoords[2 * id1 + 1];
        v2[0] = currentCoords[2 * id2];
        v2[1] = currentCoords[2 * id2 + 1];
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
    nodeCoords[0] = currentCoords[2 * nodeId];
    nodeCoords[1] = currentCoords[2 * nodeId + 1];

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
        firstNeighborCoords[0] = currentCoords[2 * firstNeighborId];
        firstNeighborCoords[1] = currentCoords[2 * firstNeighborId + 1];
        secondNeighborCoords[0] = currentCoords[2 * secondNeighborId];
        secondNeighborCoords[1] = currentCoords[2 * secondNeighborId + 1];

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
        v0[0] = currentCoords[2 * i];
        v0[1] = currentCoords[2 * i + 1];
        for (int j = 0; j < cnd; ++j) {
            int id1 = network.nodes[i].netCnxs[j];
            int id2 = network.nodes[i].netCnxs[(j + 1) % cnd];
            v1[0] = currentCoords[2 * id1];
            v1[1] = currentCoords[2 * id1 + 1];
            v2[0] = currentCoords[2 * id2];
            v2[1] = currentCoords[2 * id2 + 1];
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
            xCrds[j] = currentCoords[2 * ids[j]];
            yCrds[j] = currentCoords[2 * ids[j] + 1];
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

void LinkedNetwork::wrapCoords(std::vector<double> &coords) {
    for (int i = 0, j = 1; i < coords.size(); i += 2, j += 2) {
        coords[i] -= networkA.dimensions[0] * nearbyint(coords[i] * networkA.rpb[0]) - networkA.dimensions[0] * 0.5;
        coords[j] -= networkA.dimensions[1] * nearbyint(coords[j] * networkA.rpb[1]) - networkA.dimensions[1] * 0.5;
    }
}

// Write xyz file format of networks
void LinkedNetwork::writeXYZ(const std::string &prefix) {
    networkA.writeXYZ(prefix + "_A", "O");
    networkB.writeXYZ(prefix + "_B", "N");
}

// Write networks in format that can be loaded and visualised
void LinkedNetwork::write(const std::string &prefix) {
    networkA.write(prefix + "_A");
    networkB.write(prefix + "_B");
}
