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
    weightedDecay = inputData.weightedDecay;
    maximumBondLength = inputData.maximumBondLength;
    maximumAngle = inputData.maximumAngle * M_PI / 180; // Convert to radians
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

/**
 * @brief Populate fixedNodes with all base nodes that are a member of all fixed rings
 */
void LinkedNetwork::findFixedNodes() {
    std::for_each(fixedRings.begin(), fixedRings.end(), [this](int fixedRing) {
        for (int j = 0; j < networkB.nodes[fixedRing].dualCnxs.n; ++j) {
            if (std::find(fixedNodes.begin(), fixedNodes.end(), networkB.nodes[fixedRing].dualCnxs[j]) == fixedNodes.end())
                fixedNodes.push_back(networkB.nodes[fixedRing].dualCnxs[j]);
        }
    });
}
/**
 * @brief Calculates the difference between two vectors with periodic boundary conditions
 * @param vector1 First vector
 * @param vector2 Second vector
 * @param dimensions Dimensions of the system xhi, yhi, (zhi if 3D)
 * @throw std::invalid_argument if the sizes of the vectors do not match the dimensions
 * @return The difference vector
 */
std::vector<double> pbcVector(const std::vector<double> &vector1, const std::vector<double> &vector2, const std::vector<double> &dimensions) {
    if (vector1.size() != vector2.size() || vector1.size() != dimensions.size()) {
        std::ostringstream oss;
        oss << "Vector sizes do not match in PBC function: " << vector1.size() << " " << vector2.size() << " " << dimensions.size();
        throw std::invalid_argument(oss.str());
    }
    std::vector<double> difference_vector(vector1.size());

    for (size_t i = 0; i < vector1.size(); ++i) {
        double dimension_range = dimensions[i];
        double half_dimension_range = dimension_range / 2;
        difference_vector[i] = vector2[i] - vector1[i];
        if (difference_vector[i] > half_dimension_range) {
            difference_vector[i] -= dimension_range;
        } else if (difference_vector[i] < -half_dimension_range) {
            difference_vector[i] += dimension_range;
        }
    }
    return difference_vector;
}

/**
 * @brief read the fixed_rings.dat file and populate fixedRings with integers from each line
 * @param isFixedRingsEnabled boolean to enable or disable fixed rings
 * @param filename the name of the input file
 * @param logger the logger object
 */
void LinkedNetwork::findFixedRings(bool isFixedRingsEnabled, std::string filename, LoggerPtr logger) {
    // Format of the fixed_rings.dat file has changed to not include the number of
    // fixed rings on the first line, but simply the IDs of each fixed ring line
    // by line.
    if (isFixedRingsEnabled) {
        // Open the file
        std::ifstream fixedRingsFile(filename + ".dat", std::ios::in);
        if (!fixedRingsFile.is_open()) {
            logger->warn("Failed to open file: {}.dat, setting number of fixed rings to 0 and they will be ignored", filename);
            return;
        }
        std::string line;
        std::string ringList = "";
        // Read the file line by line
        while (getline(fixedRingsFile, line)) {
            // Convert the line to an integer and store it in the vector
            int num;
            std::istringstream(line) >> num;
            fixedRings.push_back(num);
            ringList += std::to_string(num) + " ";
        }
        logger->info("Number of fixed rings: {}", fixedRings.size());
        logger->info("Fixed rings: {}", ringList);
    } else {
        logger->info("Fixed rings disabled, setting number of fixed rings to 0.");
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

    bool foundValidMove = false;
    int baseNode1;
    int baseNode2;
    int ringNode1;
    int ringNode2;

    std::vector<int> bondBreaks;
    std::vector<int> bondMakes;
    std::vector<int> angleBreaks;
    std::vector<int> angleMakes;
    std::vector<int> ringBondBreakMake;
    std::vector<int> involvedNodes;

    for (int i = 0; i < networkA.nodes.n * networkA.nodes.n; ++i) {
        std::tuple<int, int, int, int> result;
        if (mcWeighting == "weighted")
            result = pickRandomConnection(mtGen, SelectionType::EXPONENTIAL_DECAY);
        else
            result = pickRandomConnection(mtGen, SelectionType::RANDOM);

        std::tie(baseNode1, baseNode2, ringNode1, ringNode2) = result;
        logger->debug("Picked base nodes: {} {} and ring nodes: {} {}", baseNode1, baseNode2, ringNode1, ringNode2);
        foundValidMove = genSwitchOperations(baseNode1, baseNode2, ringNode1, ringNode2,
                                             bondBreaks, bondMakes,
                                             angleBreaks, angleMakes,
                                             ringBondBreakMake, involvedNodes, logger);
        if (foundValidMove)
            break;
    }

    if (!foundValidMove) {
        logger->error("Cannot find any valid switch moves");
        throw std::runtime_error("Cannot find any valid switch moves");
    }
    numSwitches++;
    logger->info("Switch number: {}", numSwitches);

    // Save current state
    double initialEnergy = energy;

    VecF<int> saveNodeDistA = networkA.nodeDistribution;
    VecF<int> saveNodeDistB = networkB.nodeDistribution;

    VecF<VecF<int>> saveEdgeDistA = networkA.edgeDistribution;
    VecF<VecF<int>> saveEdgeDistB = networkB.edgeDistribution;
    std::vector<Node> initialInvolvedNodesA;
    for (const auto &id : involvedNodes) {
        initialInvolvedNodesA.push_back(networkA.nodes[id]);
    }
    std::vector<Node> initialInvolvedNodesB;
    for (const auto &id : ringBondBreakMake) {
        initialInvolvedNodesB.push_back(networkB.nodes[id]);
    }

    // Switch and geometry optimise
    logger->debug("Switching NetMC Network...");
    switchNetMCGraphene(bondBreaks, ringBondBreakMake);

    logger->debug("Switching LAMMPS Network...");

    if (isSimpleGrapheneEnabled) {
        logger->debug("Switching Simple Graphene");
        lammpsNetwork.switchGraphene(bondBreaks, bondMakes, angleBreaks, angleMakes, logger);
    }
    double finalEnergy;
    // Geometry optimisation of local region
    logger->debug("Minimising network...");
    lammpsNetwork.minimiseNetwork();
    std::vector<double> lammpsCoords = lammpsNetwork.getCoords(2);

    logger->info("Accepting or rejecting...");
    bool isAccepted = false;
    if (checkAnglesWithinRange(std::vector<int>{bondBreaks[0], bondBreaks[2]}, lammpsCoords, logger)) {
        if (checkBondLengths(involvedNodes, lammpsCoords, logger)) {
            finalEnergy = lammpsNetwork.getPotentialEnergy();
            if (mc.acceptanceCriterion(finalEnergy, initialEnergy, 1.0)) {
                logger->info("Accepted Move: Ei = {:.3f} Eh, Ef = {:.3f} Eh", initialEnergy, finalEnergy);
                isAccepted = true;
                numAcceptedSwitches++;
                logger->info("Syncing LAMMPS coordinates to NetMC coordinates...");
                currentCoords = lammpsCoords;
                pushCoords(lammpsCoords);
                energy = finalEnergy;
            } else {
                logger->warn("Rejected move: failed Metropolis criterion: Ei = {:.3f} Eh, Ef = {:.3f} Eh", initialEnergy, finalEnergy);
                failedEnergyChecks++;
            }
        } else {
            logger->warn("Rejected move: bond lengths are not within range");
            failedBondLengthChecks++;
        }
    } else {
        logger->warn("Rejected move: angles are not within range");
        failedAngleChecks++;
    }

    if (!isAccepted) {
        logger->debug("Reverting NetMC Network...");
        revertNetMCGraphene(initialInvolvedNodesA, initialInvolvedNodesB);

        logger->debug("Reverting LAMMPS Network...");
        lammpsNetwork.revertGraphene(bondBreaks, bondMakes, angleBreaks, angleMakes, logger);
        lammpsNetwork.setCoords(currentCoords, 2);

        networkA.nodeDistribution = saveNodeDistA;
        networkA.edgeDistribution = saveEdgeDistA;
        networkB.nodeDistribution = saveNodeDistB;
        networkB.edgeDistribution = saveEdgeDistB;
    }
}

void LinkedNetwork::showCoords(std::vector<double> &coords, LoggerPtr logger) const {
    for (int i = 0; i < coords.size(); i += 2) {
        logger->debug("{}) {} {}", i, coords[i], coords[i + 1]);
    }
}
// Rescale lattice dimensions
void LinkedNetwork::rescale(double scaleFactor) {
    networkA.rescale(scaleFactor);
    networkB.rescale(scaleFactor);
}

/**
 * @brief Chooses a random bond in the network and returns the IDs in the bond, the two rings either side and the connection type
 * @param generator Reference to a Mersenne Twister pseudo-random generator of 32-bit numbers with a state size of 19937 bits.
 * @param selectionType The type of selection to be used. Exponential decay or random.
 * @param logger The logger object.
 * @return A tuple containing the IDs of the two nodes in the bond, the two rings either side and the connection type
 * @throw std::runtime_error if the nodes in the random connection have coordinations other than 3 or 4.
 */
std::tuple<int, int, int, int> LinkedNetwork::pickRandomConnection(std::mt19937 &generator, const SelectionType &selectionType) {

    std::vector<double> weights(networkA.nodes.n);
    if (selectionType == SelectionType::EXPONENTIAL_DECAY) {
        double boxLength = dimensions[0];
        for (int i = 0; i < networkA.nodes.n; ++i) {
            double distance = networkA.nodes[i].distanceFrom(centreCoords) / boxLength;
            weights[i] = std::exp(-distance * weightedDecay);
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
        if (std::find(fixedNodes.begin(), fixedNodes.end(), randNode) == fixedNodes.end() &&
            std::find(fixedNodes.begin(), fixedNodes.end(), randNodeConnection) == fixedNodes.end()) {
            // If so, break and return the connection type
            pickingAcceptableRing = false;
        }
    }
    return std::make_tuple(randNode, randNodeConnection, sharedRingNode1, sharedRingNode2);
}

/**
 * @brief Generate the bond and angle operations for a switch move
 * @param baseNode1 the id of the first node in lattice A
 * @param baseNode2 the id of the second node in lattice A
 * @param ringNode1 the id of the first node in lattice B
 * @param ringNode2 the id of the second node in lattice B
 * @param bondBreaks the bonds to break
 * @param bondMakes the bonds to make
 * @param angleBreaks the angles to break
 * @param angleMakes the angles to make
 * @param ringBondBreakMake the ring bond to break and make, first two indexes are bond to break, second two are bond to make
 * @return true if the switch move is possible, false otherwise
 */
bool LinkedNetwork::genSwitchOperations(int baseNode1, int baseNode2, int ringNode1, int ringNode2,
                                        std::vector<int> &bondBreaks, std::vector<int> &bondMakes,
                                        std::vector<int> &angleBreaks, std::vector<int> &angleMakes,
                                        std::vector<int> &ringBondBreakMake, std::vector<int> &convexCheckIDs, LoggerPtr logger) {
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
     *               \       /        (6 membered case)     |      |
     *                9-----10                              9------10
     *
     *      Bonds to break       Bonds to Make
     *      1-5, 2-4             1-4, 2-5
     *
     *      Angles to break      Angles to Make
     *      1-5-9, 1-5-13        1-4-8, 1-4-12
     *      2-4-8, 2-4-12        2-5-9, 2-5-13
     *      4-2-1, 4-2-6         4-1-2, 4-1-3
     *      5-1-2, 5-1-3         1-2-5, 6-2-5
     */
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

    // Additional error checking
    if (baseNode5 == baseNode6 || baseNode3 == baseNode4) {
        logger->debug("No valid move: Selected nodes describe an edge of two edge sharing triangles");
        return false;
    }
    // Prevent rings having only two or fewer neighbours
    if (getUniqueValues(networkB.nodes[ringNode1].netCnxs).n <= 3 || getUniqueValues(networkB.nodes[ringNode2].netCnxs).n <= 3) {
        logger->debug("No valid move: Switch would result in a ring size less than 3");
        return false;
    }

    // check move will not violate dual connectivity limits
    if (networkB.nodes[ringNode1].netCnxs.n == minBCnxs ||
        networkB.nodes[ringNode2].netCnxs.n == minBCnxs ||
        networkB.nodes[ringNode3].netCnxs.n == maxBCnxs ||
        networkB.nodes[ringNode4].netCnxs.n == maxBCnxs) {
        logger->debug("No valid move: Switch would violate dual connectivity limits");
        return false;
    }
    bondBreaks = {baseNode1, baseNode5, baseNode2, baseNode4};

    bondMakes = {baseNode1, baseNode4, baseNode2, baseNode5};
    //                          Break                   Make
    ringBondBreakMake = {ringNode1, ringNode2, ringNode3, ringNode4};

    convexCheckIDs = {baseNode1, baseNode2, baseNode3, baseNode4, baseNode5,
                      baseNode6, baseNode7, baseNode8, baseNode9, baseNode10,
                      baseNode11, baseNode12, baseNode13, baseNode14};

    if (baseNode7 == baseNode4) {
        // 4 membered ringNode2
        angleBreaks = {baseNode3, baseNode4, baseNode2,
                       baseNode12, baseNode4, baseNode2,
                       baseNode1, baseNode2, baseNode4,
                       baseNode6, baseNode2, baseNode4};
        angleMakes = {baseNode3, baseNode4, baseNode1,
                      baseNode12, baseNode4, baseNode1,
                      baseNode3, baseNode1, baseNode4,
                      baseNode2, baseNode1, baseNode4};
        logger->debug(" {:03}------{:03}------{:03}------{:03}             {:03}-----{:03}-----{:03}-----{:03} ", baseNode11, baseNode3, baseNode4, baseNode12, baseNode11, baseNode3, baseNode4, baseNode12);
        logger->debug("            |       |                                \\ {:03}  / ", ringNode2, ringNode2);
        logger->debug("            |   {:03}  |                             \\    /");
        logger->debug("            |       |                                  {:03}", baseNode1);

    } else if (baseNode7 == baseNode8) {
        // 5 membered ringNode2
        angleBreaks = {baseNode7, baseNode4, baseNode2,
                       baseNode12, baseNode4, baseNode2,
                       baseNode1, baseNode2, baseNode4,
                       baseNode6, baseNode2, baseNode4};
        angleMakes = {baseNode7, baseNode4, baseNode1,
                      baseNode12, baseNode4, baseNode1,
                      baseNode2, baseNode1, baseNode4,
                      baseNode3, baseNode1, baseNode4};
        logger->debug("");
        logger->debug("               {:03}                                      {:03}", baseNode7, baseNode7);
        logger->debug("            /      \\                                   /   \\");
        logger->debug("           /        \\                                 /     \\");
        logger->debug(" {:03}-----{:03}   {:03}  {:03}-----{:03}             {:03}-----{:03} {:03} {:03}-----{:03}", baseNode11, baseNode3, ringNode2, baseNode4, baseNode12, baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
        logger->debug("          \\          /                                \\     /");
        logger->debug("           \\        /                                   {:03}", baseNode1);

    } else {
        // 6+ membered ringNode2
        angleBreaks = {baseNode2, baseNode4, baseNode8,
                       baseNode2, baseNode4, baseNode12,
                       baseNode4, baseNode2, baseNode1,
                       baseNode4, baseNode2, baseNode6};
        angleMakes = {baseNode1, baseNode4, baseNode8,
                      baseNode1, baseNode4, baseNode12,
                      baseNode4, baseNode1, baseNode2,
                      baseNode4, baseNode1, baseNode3};
        logger->debug("");
        logger->debug("           {:03}~~~~~{:03}                              {:03}~~~~~{:03}", baseNode7, baseNode8, baseNode7, baseNode8);
        logger->debug("           /        \\                                |       |");
        logger->debug("          /          \\                      {:03}-----{:03} {:03} {:03}-----{:03}", baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
        logger->debug(" {:03}-----{:03}   {:03}   {:03}-----{:03}                      \\     /", baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
        logger->debug("          \\          /                                 \\   /");
        logger->debug("           \\        /                                   {:03}", baseNode1);
    }
    logger->debug("    {:03}    {:03}-----{:03}   {:03}          --->     {:03}       |      {:03}", ringNode3, baseNode1, baseNode2, ringNode4, ringNode3, ringNode4);
    if (baseNode5 == baseNode10) {
        // 4 membered ringNode1
        angleBreaks.insert(angleBreaks.end(), {baseNode13, baseNode5, baseNode1,
                                               baseNode6, baseNode5, baseNode1,
                                               baseNode3, baseNode1, baseNode5,
                                               baseNode2, baseNode1, baseNode5});
        angleMakes.insert(angleMakes.end(), {baseNode13, baseNode5, baseNode2,
                                             baseNode6, baseNode5, baseNode2,
                                             baseNode1, baseNode2, baseNode5,
                                             baseNode6, baseNode2, baseNode5});
        logger->debug("            |       |                                   {:03}", baseNode2);
        logger->debug("            |       |                                  /   \\");
        logger->debug("            |  {:03}  |                                 / {:03} \\ ", ringNode1, ringNode1);
        logger->debug(" {:03}-------{:03}-----{:03}-------{:03}            {:03}-----{:03}-----{:03}-----{:03} ", baseNode13, baseNode5, baseNode6, baseNode14, baseNode13, baseNode5, baseNode6, baseNode14);
        logger->debug("");
    } else if (baseNode9 == baseNode10) {
        // 5 membered ringNode1
        angleBreaks.insert(angleBreaks.end(), {baseNode13, baseNode5, baseNode1,
                                               baseNode9, baseNode5, baseNode1,
                                               baseNode3, baseNode1, baseNode5,
                                               baseNode2, baseNode1, baseNode5});
        angleMakes.insert(angleMakes.end(), {baseNode13, baseNode5, baseNode2,
                                             baseNode9, baseNode5, baseNode2,
                                             baseNode1, baseNode2, baseNode5,
                                             baseNode6, baseNode2, baseNode5});

        logger->debug("           /        \\                                   {:03}", baseNode2);
        logger->debug("          /          \\                                /     \\");
        logger->debug(" {:03}-----{:03}   {:03}  {:03}-----{:03}             {:03}-----{:03} {:03} {:03}-----{:03}", baseNode13, baseNode5, ringNode1, baseNode6, baseNode14, baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
        logger->debug("           \\        /                                 \\     /");
        logger->debug("            \\      /                                   \\   /");
        logger->debug("              {:03}                                       {:03}", baseNode9, baseNode9);
        logger->debug("");
    } else {
        // 6+ membered ringNode1
        angleBreaks.insert(angleBreaks.end(), {baseNode1, baseNode5, baseNode9,
                                               baseNode1, baseNode5, baseNode13,
                                               baseNode5, baseNode1, baseNode2,
                                               baseNode5, baseNode1, baseNode3});
        angleMakes.insert(angleMakes.end(), {baseNode2, baseNode5, baseNode9,
                                             baseNode2, baseNode5, baseNode13,
                                             baseNode1, baseNode2, baseNode5,
                                             baseNode6, baseNode2, baseNode5});
        logger->debug("           /        \\                                   {:03}", baseNode2);
        logger->debug("          /          \\                                 /   \\");
        logger->debug(" {:03}-----{:03}   {:03}   {:03}-----{:03}                      /     \\", baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
        logger->debug("          \\          /                      {:03}-----{:03} {:03} {:03}-----{:03}", baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
        logger->debug("           \\        /                                |       |");
        logger->debug("           {:03}~~~~~{:03}                              {:03}~~~~~{:03}", baseNode9, baseNode10, baseNode9, baseNode10);
        logger->debug("");
    }

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
    VecR<int> commonConnections = vCommonValues(networkA.nodes[baseNode].netCnxs, networkB.nodes[ringNode].dualCnxs);
    commonConnections.delValue(excludeNode);
    if (commonConnections.n == 1)
        return commonConnections[0];

    std::ostringstream oss;
    oss << "Could not find common base node for base node " << baseNode << " and ring node " << ringNode
        << " excluding base node " << excludeNode << " common connections: " << commonConnections;
    throw std::runtime_error(oss.str());
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
void LinkedNetwork::switchNetMCGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &ringBondBreakMake) {
    if (bondBreaks.size() != 4 || ringBondBreakMake.size() != 4) {
        throw std::invalid_argument("Invalid input sizes for switchNetMCGraphene");
    }

    int atom1 = bondBreaks[0];
    int atom2 = bondBreaks[2];
    int atom4 = bondBreaks[3];
    int atom5 = bondBreaks[1];

    int ringNode1 = ringBondBreakMake[0];
    int ringNode2 = ringBondBreakMake[1];
    int ringNode3 = ringBondBreakMake[2];
    int ringNode4 = ringBondBreakMake[3];

    // Apply changes to descriptors due to breaking connections
    // For network A node distribution and edge distribution will remain unchanged

    // logger->debug("Switching A-A connections");
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
    networkB.nodes[ringNode3].netCnxs.addValue(ringNode4);
    networkB.nodes[ringNode4].netCnxs.addValue(ringNode3);

    // logger->debug("Switching B-A connections");
    // B-A connectivities
    networkB.nodes[ringNode1].dualCnxs.delValue(atom1);
    networkB.nodes[ringNode2].dualCnxs.delValue(atom2);
    networkB.nodes[ringNode3].dualCnxs.addValue(atom2);
    networkB.nodes[ringNode4].dualCnxs.addValue(atom1);
}

void LinkedNetwork::revertNetMCGraphene(const std::vector<Node> &initialInvolvedNodesA, const std::vector<Node> &initialInvolvedNodesB) {
    // Revert changes to descriptors due to breaking connections
    for (const Node &node : initialInvolvedNodesA) {
        networkA.nodes[node.id] = node;
    }
    for (const Node &node : initialInvolvedNodesB) {
        networkB.nodes[node.id] = node;
    }
}

/**
 * @brief Checks if a ring's base nodes have more than two rings in common
 * @param ringNode ID of the ring node
 * @return true if the ring's base nodes have at most two rings in common, false otherwise
 */
bool LinkedNetwork::checkThreeRingEdges(const int &ringNode) const {
    int ringSize = networkB.nodes[ringNode].dualCnxs.n;
    for (int i = 0; i < ringSize; ++i) {
        int baseNode1 = networkB.nodes[ringNode].dualCnxs[i];
        int baseNode2 = networkB.nodes[ringNode].dualCnxs[(i + 1) % ringSize];
        if (vCommonValues(networkA.nodes[baseNode1].dualCnxs, networkA.nodes[baseNode2].dualCnxs).n > 2) {
            // The two base nodes have two or more rings in common, this is not possible unless they are the same node.
            return false;
        }
    }
    return true;
}

/**
 * @brief Sets the coordinates of nodes in the network. Ring network coordinates are recalculated.
 * @param coords Coordinates of the base nodes
 * @throw std::invalid_argument if the number of coordinates does not match the number of nodes in network A
 */
void LinkedNetwork::pushCoords(std::vector<double> &coords) {
    if (coords.size() != 2 * networkA.nodes.n) {
        throw std::invalid_argument("Number of coordinates does not match number of nodes in network A");
    }

    // Sync A coordinates
    for (int i = 0; i < networkA.nodes.n; ++i) {
        networkA.nodes[i].crd[0] = coords[2 * i];
        networkA.nodes[i].crd[1] = coords[2 * i + 1];
    }
    // Sync B coordinates
    std::vector<double> total(2, 0.0);
    for (int ringNode = 0; ringNode < networkB.nodes.n; ++ringNode) {
        Node &selectedRingNode = networkB.nodes[ringNode];
        total[0] = total[1] = 0.0;
        for (int neighbour = 0; neighbour < selectedRingNode.dualCnxs.n; ++neighbour) {
            std::vector<double> pbcCoords = pbcVector(selectedRingNode.crd, networkA.nodes[selectedRingNode.dualCnxs[neighbour]].crd, dimensions);
            total[0] += pbcCoords[0];
            total[1] += pbcCoords[1];
        }
        total[0] /= selectedRingNode.dualCnxs.n;
        total[1] /= selectedRingNode.dualCnxs.n;

        selectedRingNode.crd[0] += total[0];
        selectedRingNode.crd[1] += total[1];

        // Wrap the new coordinates back into the box
        for (size_t i = 0; i < selectedRingNode.crd.n; ++i) {
            while (selectedRingNode.crd[i] < 0) {
                selectedRingNode.crd[i] += dimensions[i];
            }
            while (selectedRingNode.crd[i] >= dimensions[i]) {
                selectedRingNode.crd[i] -= dimensions[i];
            }
        }
    }
}
// Get normalised probability distribution of nodes of each size in given
// lattice
VecF<double> LinkedNetwork::getNodeDistribution(std::string_view networkType) {
    if (networkType == "A")
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
VecF<int> LinkedNetwork::getMaxClusters(std::string_view lattice, int minCnd, int maxCnd) {
    if (lattice == "A")
        return networkA.maxClusters(minCnd, maxCnd, 3, 2);
    else
        return networkB.maxClusters(minCnd, maxCnd, 3, 2);
}

// Check linked networks for consistency
bool LinkedNetwork::checkConsistency() {
    return (checkCnxConsistency() && checkDescriptorConsistency());
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

// Get sum of areas and squared areas for each ring size
void LinkedNetwork::getRingAreas(VecF<double> &areaSum, VecF<double> &areaSqSum) {

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
/**
 * @brief Wraps the coordinates out of bounds back into the periodic box, only for 2 dimensions
 * @param coords Coordinates to be wrapped (1D vector of pairs)
 */
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

/**
 * @brief Calculates angle between two vectors
 * @param vector1 First vector
 * @param vector2 Second vector
 * @return Angle between the two vectors in radians
 */
double getClockwiseAngleBetweenVectors(const std::vector<double> &vector1, const std::vector<double> &vector2) {
    double dotProduct = vector1[0] * vector2[0] + vector1[1] * vector2[1];
    double magnitudeProduct = std::sqrt(vector1[0] * vector1[0] + vector1[1] * vector1[1]) *
                              std::sqrt(vector2[0] * vector2[0] + vector2[1] * vector2[1]);
    double angle = std::acos(dotProduct / magnitudeProduct);

    // Calculate the cross product of the two vectors
    double crossProduct = vector1[0] * vector2[1] - vector1[1] * vector2[0];

    // If the cross product is positive, subtract the angle from 2π to get the angle in the range [π, 2π]
    if (crossProduct > 0) {
        angle = 2 * M_PI - angle;
    }

    return angle;
}

/**
 *
 * @brief Get the clockwise angle of the vector between two nodes relative to the x axis
 * @param coord1 Coordinates of the first node
 * @param coord2 Coordinates of the second node
 * @param dimensions Dimensions of the periodic box
 * @return Angle of the vector between the two nodes relative to the x axis in radians
 */
double getClockwiseAngle(const std::vector<double> &coord1, const std::vector<double> &coord2, const std::vector<double> &dimensions) {
    std::vector<double> periodicVector = pbcVector(coord1, coord2, dimensions);

    // Range of this function is -pi to pi
    double angle = std::atan2(periodicVector[1], periodicVector[0]);

    // So we have to convert it to 0 to 2pi
    if (angle < 0) {
        angle += 2 * M_PI;
    }
    return 2 * M_PI - angle;
}

/**
 * @brief Checks if a given node has clockwise neighbours
 * @param nodeID ID of the node to check
 * @return true if the node has clockwise neighbours, false otherwise
 */
bool LinkedNetwork::checkClockwiseNeighbours(const int &nodeID) const {
    Node node = networkA.nodes[nodeID];
    double prevAngle = getClockwiseAngle(node.crd, networkA.nodes[node.netCnxs.back()].crd, dimensions);

    // The angle a neighbour makes with the x axis should always increase, EXCEPT when the angle wraps around from 2π to 0
    // This should only happen once if the neighbours are clockwise
    int timesDecreased = 0;
    for (int id : node.netCnxs) {
        double angle = getClockwiseAngle(node.crd, networkA.nodes[id].crd, dimensions);
        if (angle < prevAngle) {
            timesDecreased++;
            if (timesDecreased == 2) {
                // So return false if it decreases twice
                return false;
            }
        }
        prevAngle = angle;
    }
    return true;
}
/**
 * @brief Checks if a given node in network A has clockwise neighbours with given coordinates
 * @param nodeID ID of the node to check
 * @return true if the node has clockwise neighbours, false otherwise
 */
bool LinkedNetwork::checkClockwiseNeighbours(const int &nodeID, const std::vector<double> &coords) const {
    const std::vector<double> nodeCoord = {coords[2 * nodeID], coords[2 * nodeID + 1]};
    const int lastNeighbourCoordsID = networkA.nodes[nodeID].netCnxs.back();
    const std::vector<double> lastNeighbourCoord = {coords[2 * lastNeighbourCoordsID], coords[2 * lastNeighbourCoordsID + 1]};
    double prevAngle = getClockwiseAngle(nodeCoord, lastNeighbourCoord, dimensions);
    int timesDecreased = 0;
    for (const auto &id : networkA.nodes[nodeID].netCnxs) {
        double angle = getClockwiseAngle(nodeCoord, {coords[2 * id], coords[2 * id + 1]}, dimensions);
        if (angle < prevAngle) {
            timesDecreased++;
            if (timesDecreased == 2) {
                return false;
            }
        }
        prevAngle = angle;
    }
    return true;
}

/**
 * @brief Checks if all nodes have clockwise neighbours
 * @param logger Logger to log any anticlockwise neighbours
 * @return true if all nodes have clockwise neighbours, false otherwise
 */
bool LinkedNetwork::checkAllClockwiseNeighbours(LoggerPtr logger) const {
    bool check = true;
    for (int nodeID = 0; nodeID < networkA.nodes.n; ++nodeID) {
        if (!checkClockwiseNeighbours(nodeID)) {
            std::ostringstream oss;
            oss << "Node " << nodeID << " has anticlockwise neighbours: ";
            for (int neighbourID : networkA.nodes[nodeID].netCnxs) {
                oss << neighbourID << " ";
            }
            logger->warn(oss.str());
            check = false;
        }
    }
    return check;
}

void LinkedNetwork::arrangeNeighboursClockwise(const int &nodeID, const std::vector<double> &coords) {
    // Get the coordinates of the center node
    const std::vector<double> nodeCoord = {coords[2 * nodeID], coords[2 * nodeID + 1]};

    // Get the neighbour IDs of the center node
    VecR<int> &neighbours = networkA.nodes[nodeID].netCnxs;

    // Create a vector of pairs, where each pair contains a neighbour ID and the angle between the center node and that neighbour
    std::vector<std::pair<int, double>> neighbourAngles;
    for (const auto &neighbourID : neighbours) {
        std::vector<double> neighbourCoord = {coords[2 * neighbourID], coords[2 * neighbourID + 1]};
        double angle = getClockwiseAngle(nodeCoord, neighbourCoord, dimensions);
        neighbourAngles.push_back({neighbourID, angle});
    }

    // Sort the vector of pairs in ascending order of angle
    std::sort(neighbourAngles.begin(), neighbourAngles.end(), [](const std::pair<int, double> &a, const std::pair<int, double> &b) {
        return a.second < b.second;
    });

    // Replace the neighbour IDs in the center node with the sorted neighbour IDs
    for (size_t i = 0; i < neighbours.n; ++i) {
        neighbours[i] = neighbourAngles[i].first;
    }
}

/**
 * @brief Checks if all angles around all nodes are less than or equal to maximum angle
 * @param coords Coordinates of all nodes as a 1D vector of coordinate pairs
 * @param logger Logger to log any non-convex angles
 * @return true if all angles are convex, false otherwise
 */
bool LinkedNetwork::checkAnglesWithinRange(const std::vector<double> &coords, LoggerPtr logger) {
    for (int nodeID = 0; nodeID < networkA.nodes.n; ++nodeID) {
        arrangeNeighboursClockwise(nodeID, coords);
        for (int i = 0; i < networkA.nodes[nodeID].netCnxs.n; ++i) {
            int neighbourID = networkA.nodes[nodeID].netCnxs[i];
            int nextNeighbourID = networkA.nodes[nodeID].netCnxs[(i + 1) % networkA.nodes[nodeID].netCnxs.n];
            std::vector<double> v1 = pbcVector(std::vector<double>{coords[neighbourID * 2], coords[neighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            std::vector<double> v2 = pbcVector(std::vector<double>{coords[nextNeighbourID * 2], coords[nextNeighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            double angle = getClockwiseAngleBetweenVectors(v1, v2);
            if (angle > maximumAngle) {
                logger->warn("Node {} has an out of bounds angle between neighbours {} and {}: {:.2f} degrees", nodeID, neighbourID, nextNeighbourID, angle * 180 / M_PI);

                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Checks if all angles around the given nodes are less than or equal to maximum angle
 * @param nodeIDs IDs of the nodes to check
 * @param coords Coordinates of all nodes as a 1D vector of coordinate pairs
 * @param logger Logger to log any non-convex angles
 * @return true if all angles are convex, false otherwise
 */
bool LinkedNetwork::checkAnglesWithinRange(const std::vector<int> &nodeIDs, const std::vector<double> &coords, LoggerPtr logger) {
    for (int nodeID : nodeIDs) {
        // if (std::find(fixedNodes.begin(), fixedNodes.end(), nodeID) != fixedNodes.end()) {
        arrangeNeighboursClockwise(nodeID, coords);
        for (int i = 0; i < networkA.nodes[nodeID].netCnxs.n; ++i) {
            int neighbourID = networkA.nodes[nodeID].netCnxs[i];
            int nextNeighbourID = networkA.nodes[nodeID].netCnxs[(i + 1) % networkA.nodes[nodeID].netCnxs.n];
            std::vector<double> v1 = pbcVector(std::vector<double>{coords[neighbourID * 2], coords[neighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            std::vector<double> v2 = pbcVector(std::vector<double>{coords[nextNeighbourID * 2], coords[nextNeighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            double angle = getClockwiseAngleBetweenVectors(v1, v2);
            if (angle > maximumAngle) {
                logger->warn("Node {} has an out of bounds angle between neighbours {} and {}: {:.2f} degrees", nodeID, neighbourID, nextNeighbourID, angle * 180 / M_PI);
                return false;
            }
        }
    }
    //}
    return true;
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

bool LinkedNetwork::checkBondLengths(const int &nodeID, const std::vector<double> &coords, LoggerPtr logger) const {
    for (const auto &neighbourID : networkA.nodes[nodeID].netCnxs) {
        std::vector<double> pbcVec = pbcVector(std::vector<double>{coords[2 * nodeID], coords[2 * nodeID + 1]},
                                               std::vector<double>{coords[2 * neighbourID], coords[2 * neighbourID + 1]}, dimensions);
        if (std::sqrt(pbcVec[0] * pbcVec[0] + pbcVec[1] * pbcVec[1]) > maximumBondLength) {
            logger->warn("Node {} has a bond length greater than the maximum bond length with neighbour {}: {:.2f}", nodeID, neighbourID, std::sqrt(pbcVec[0] * pbcVec[0] + pbcVec[1] * pbcVec[1]));
            return false;
        }
    }
    return true;
}

bool LinkedNetwork::checkBondLengths(const std::vector<int> &nodeIDs, const std::vector<double> &coords, LoggerPtr logger) const {
    return std::all_of(nodeIDs.begin(), nodeIDs.end(), [this, &coords, &logger](int nodeID) {
        return checkBondLengths(nodeID, coords, logger);
    });
}