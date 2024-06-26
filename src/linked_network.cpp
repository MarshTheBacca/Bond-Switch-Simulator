#include "linked_network.h"
#include <filesystem>

/**
 * @brief Default constructor
 */
LinkedNetwork::LinkedNetwork() = default;


/**
 * @brief Construct by loading networks from files
 * @param inputData the input data object
 * @param loggerArg the logger object
 */
LinkedNetwork::LinkedNetwork(const InputData &inputData, const LoggerPtr &loggerArg) : minRingSize(inputData.minRingSize),
                                                                                       maxRingSize(inputData.maxRingSize),
                                                                                       selectionType(inputData.randomOrWeighted),
                                                                                       metropolisCondition(inputData.randomSeed),
                                                                                       weightedDecay(inputData.weightedDecay),
                                                                                       maximumBondLength(inputData.maximumBondLength),
                                                                                       maximumAngle(inputData.maximumAngle * M_PI / 180),
                                                                                       writeMovie(inputData.writeMovie),
                                                                                       logger(loggerArg) {
    networkA = Network(NetworkType::BASE_NETWORK, logger);
    networkB = Network(NetworkType::DUAL_NETWORK, logger);

    if (inputData.isFixRingsEnabled) {
        findFixedRings(std::filesystem::path("./input_files") / "bss_network" /"fixed_rings.txt");
        findFixedNodes();
    } else {
        logger->info("Fixed rings disabled, setting number of fixed rings to 0.");
    }
    if (int loadedMinRingSize = networkB.getMinConnections(); loadedMinRingSize < minRingSize) {
        logger->warn("Loaded network has a min ring size of {} which is lower than input file's {}", loadedMinRingSize, minRingSize);
    }
    if (int loadedMaxRingSize = networkB.getMaxConnections(); loadedMaxRingSize > maxRingSize) {
        logger->warn("Loaded network has a max ring size of {} which is higher than input file's {}", loadedMaxRingSize, maxRingSize);
    }
    dimensions = networkA.dimensions;
    centreCoords = {dimensions[0] / 2, dimensions[1] / 2};

    lammpsNetwork = LammpsObject(logger);
    if (writeMovie) {
        lammpsNetwork.startMovie();
        lammpsNetwork.writeMovie();
    }
    lammpsNetwork.minimiseNetwork();
    currentCoords = lammpsNetwork.getCoords(2);
    energy = lammpsNetwork.getPotentialEnergy();
    pushCoords(currentCoords);
    weights.resize(networkA.nodes.size());
    updateWeights();
    randomNumGen.seed(inputData.randomSeed);
}

/**
 * @brief read the fixed_rings.txt file and populate fixedRings with integers from each line
 * @param isFixedRingsEnabled boolean to enable or disable fixed rings
 * @param filename the name of the input file
 */
void LinkedNetwork::findFixedRings(const std::string &filePath) {
    // File structure has changed to have each fixed ring per line
    // ie, file does not start with the number of fixed rings.
    std::ifstream fixedRingsFile(filePath, std::ios::in);
    if (!fixedRingsFile.is_open()) {
        logger->warn("Failed to open file: {}, setting number of fixed rings to 0 and they will be ignored", filePath);
        return;
    }
    std::string line;
    while (std::getline(fixedRingsFile, line)) {
        int num;
        std::istringstream(line) >> num;
        fixedRings[num] = networkB.nodes[num].netConnections.size();
    }
    logger->info("Number of fixed rings: {}", fixedRings.size());
    logger->info("Fixed ring info (ID: Size): {}", mapToString(fixedRings));
}

/**
 * @brief Populate fixedNodes with all base nodes that are a member of all fixed rings
 */
void LinkedNetwork::findFixedNodes() {
    std::for_each(fixedRings.begin(), fixedRings.end(), [this](const std::pair<const int, int> &fixedRingPair) {
        int fixedRing = fixedRingPair.first;
        std::for_each(networkB.nodes[fixedRing].dualConnections.begin(), networkB.nodes[fixedRing].dualConnections.end(), [this](int fixedNode) {
            fixedNodes.insert(fixedNode);
        });
    });
}

/**
 * @brief Perform a monte carlo switch move, evaluate energy, and accept or reject
 */
void LinkedNetwork::monteCarloSwitchMoveLAMMPS(const double &temperature) {
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
    std::unordered_set<int> involvedNodes;

    for (int i = 0; i < networkA.nodes.size() * networkA.nodes.size(); ++i) {
        std::tuple<int, int, int, int> result;
        result = pickRandomConnection();

        std::tie(baseNode1, baseNode2, ringNode1, ringNode2) = result;
        logger->debug("Picked base nodes: {} {} and ring nodes: {} {}", baseNode1, baseNode2, ringNode1, ringNode2);
        foundValidMove = genSwitchOperations(baseNode1, baseNode2, ringNode1, ringNode2,
                                             bondBreaks, bondMakes,
                                             angleBreaks, angleMakes,
                                             ringBondBreakMake, involvedNodes);
        if (foundValidMove)
            break;
    }

    if (!foundValidMove) {
        logger->error("Cannot find any valid switch moves");
        throw std::runtime_error("Cannot find any valid switch moves");
    }
    numSwitches++;
    logger->debug("Switch number: {}", numSwitches);

    // Save current state
    double initialEnergy = energy;

    std::vector<Node> initialInvolvedNodesA;
    for (const auto &id : involvedNodes) {
        initialInvolvedNodesA.push_back(networkA.nodes[id]);
    }
    std::vector<Node> initialInvolvedNodesB;
    for (const auto &id : ringBondBreakMake) {
        initialInvolvedNodesB.push_back(networkB.nodes[id]);
    }

    // Switch and geometry optimise
    logger->debug("Switching BSS Network...");
    switchNetMCGraphene(bondBreaks, ringBondBreakMake);
    std::vector<double> rotatedCoord1;
    std::vector<double> rotatedCoord2;
    std::vector<int> orderedRingNodes = {ringBondBreakMake[1], ringBondBreakMake[3], ringBondBreakMake[0], ringBondBreakMake[2]};
    std::tie(rotatedCoord1, rotatedCoord2) = rotateBond(baseNode1, baseNode2, getRingsDirection(orderedRingNodes));

    logger->debug("Switching LAMMPS Network...");

    lammpsNetwork.switchGraphene(bondBreaks, bondMakes, angleBreaks, angleMakes, rotatedCoord1, rotatedCoord2);

    // Geometry optimisation of local region
    logger->debug("Minimising network...");
    lammpsNetwork.minimiseNetwork();
    std::vector<double> lammpsCoords = lammpsNetwork.getCoords(2);

    logger->debug("Accepting or rejecting...");
    if (!checkAnglesWithinRange(setDifference(involvedNodes, fixedNodes), lammpsCoords)) {
        logger->debug("Rejected move: angles are not within range");
        failedAngleChecks++;
        rejectMove(initialInvolvedNodesA, initialInvolvedNodesB, bondBreaks, bondMakes, angleBreaks, angleMakes);
        return;
    }
    if (!checkBondLengths(involvedNodes, lammpsCoords)) {
        logger->debug("Rejected move: bond lengths are not within range");
        failedBondLengthChecks++;
        rejectMove(initialInvolvedNodesA, initialInvolvedNodesB, bondBreaks, bondMakes, angleBreaks, angleMakes);
        return;
    }
    double finalEnergy = lammpsNetwork.getPotentialEnergy();
    if (!metropolisCondition.acceptanceCriterion(finalEnergy, initialEnergy, temperature)) {
        logger->debug("Rejected move: failed Metropolis criterion: Ei = {:.3f} Eh, Ef = {:.3f} Eh", initialEnergy, finalEnergy);
        failedEnergyChecks++;
        rejectMove(initialInvolvedNodesA, initialInvolvedNodesB, bondBreaks, bondMakes, angleBreaks, angleMakes);
        return;
    }
    logger->debug("Accepted Move: Ei = {:.3f} Eh, Ef = {:.3f} Eh", initialEnergy, finalEnergy);
    numAcceptedSwitches++;
    logger->debug("Syncing LAMMPS coordinates to BSS coordinates...");
    currentCoords = lammpsCoords;
    pushCoords(currentCoords);
    updateWeights();
    arrangeNeighboursClockwise(involvedNodes, currentCoords);
    energy = finalEnergy;
    if (writeMovie)
        lammpsNetwork.writeMovie();
}

void LinkedNetwork::rejectMove(const std::vector<Node> &initialInvolvedNodesA, const std::vector<Node> &initialInvolvedNodesB,
                               const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                               const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes) {
    logger->debug("Reverting BSS Network...");
    revertNetMCGraphene(initialInvolvedNodesA, initialInvolvedNodesB);
    logger->debug("Reverting LAMMPS Network...");
    lammpsNetwork.revertGraphene(bondBreaks, bondMakes, angleBreaks, angleMakes);
    logger->debug("Syncing LAMMPS coordinates to BSS coordinates...");
    lammpsNetwork.setCoords(currentCoords, 2);
}

void LinkedNetwork::showCoords(const std::vector<double> &coords) const {
    for (int i = 0; i < coords.size(); i += 2) {
        logger->debug("{}) {} {}", i, coords[i], coords[i + 1]);
    }
}

/**
 * @brief Scale the entire linked network by a given factor
 * @param scaleFactor the factor to scale the network by
 */
void LinkedNetwork::rescale(double scaleFactor) {
    vectorMultiply(dimensions, scaleFactor);
    networkA.rescale(scaleFactor);
    networkB.rescale(scaleFactor);
}

void LinkedNetwork::updateWeights() {
    if (selectionType == SelectionType::EXPONENTIAL_DECAY) {
        double boxLength = dimensions[0];
        for (int i = 0; i < networkA.nodes.size(); ++i) {
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
}

/**
 * @brief Chooses a random bond in the network and returns the IDs in the bond, the two rings either side and the connection type
 * @return A tuple containing the IDs of the two nodes in the bond, the two rings either side and the connection type
 * @throw std::runtime_error if the nodes in the random connection have coordinations other than 3 or 4.
 */
std::tuple<int, int, int, int> LinkedNetwork::pickRandomConnection() {

    // Create a discrete distribution based on the weights
    std::discrete_distribution<> distribution(weights.begin(), weights.end());
    int randNode;
    int randNodeConnection;
    int sharedRingNode1;
    int sharedRingNode2;
    bool pickingAcceptableRing = true;
    std::uniform_int_distribution<int> randomCnx;
    std::uniform_int_distribution randomDirection(0, 1);

    while (pickingAcceptableRing) {
        randNode = distribution(randomNumGen);
        int randNodeCoordination = networkA.nodes[randNode].netConnections.size();
        randomCnx.param(std::uniform_int_distribution<int>::param_type(0, randNodeCoordination - 1));
        randNodeConnection = networkA.nodes[randNode].netConnections[randomCnx(randomNumGen)];

        // Two connected nodes should always share two ring nodes.
        if (std::unordered_set<int> commonRings = intersectVectors(networkA.nodes[randNode].dualConnections,
                                                                   networkA.nodes[randNodeConnection].dualConnections);
            commonRings.size() == 2) {
            // Randomly assign ringNode1 and ringNode2 to those two common ring nodes
            auto it = commonRings.begin();
            std::advance(it, randomDirection(randomNumGen));
            sharedRingNode1 = *it;
            commonRings.erase(it);
            sharedRingNode2 = *commonRings.begin();
        } else {
            logger->warn("Selected random connection does not share two ring nodes: {} {}", randNode, randNodeConnection);
            continue;
        }
        pickingAcceptableRing = false;
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
                                        std::vector<int> &ringBondBreakMake, std::unordered_set<int> &convexCheckIDs) {
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
    int baseNode5 = findCommonConnection(baseNode1, ringNode1, baseNode2);
    int baseNode6 = findCommonConnection(baseNode2, ringNode1, baseNode1);
    int baseNode3 = findCommonConnection(baseNode1, ringNode2, baseNode2);
    int baseNode4 = findCommonConnection(baseNode2, ringNode2, baseNode1);

    int ringNode3 = findCommonRing(baseNode1, baseNode5, ringNode1);
    int ringNode4 = findCommonRing(baseNode2, baseNode6, ringNode1);

    int baseNode11 = findCommonConnection(baseNode3, ringNode3, baseNode1);
    int baseNode7 = findCommonConnection(baseNode3, ringNode2, baseNode1);
    int baseNode8 = findCommonConnection(baseNode4, ringNode2, baseNode2);
    int baseNode12 = findCommonConnection(baseNode4, ringNode4, baseNode2);
    int baseNode14 = findCommonConnection(baseNode6, ringNode4, baseNode2);
    int baseNode10 = findCommonConnection(baseNode6, ringNode1, baseNode2);
    int baseNode9 = findCommonConnection(baseNode5, ringNode1, baseNode1);
    int baseNode13 = findCommonConnection(baseNode5, ringNode3, baseNode1);

    // Additional error checking
    if (baseNode5 == baseNode6 || baseNode3 == baseNode4) {
        logger->debug("No valid move: Selected nodes describe an edge of two edge sharing triangles");
        return false;
    }
    // Prevent rings having only two or fewer neighbours
    if (networkB.nodes[ringNode1].netConnections.size() <= 3 || networkB.nodes[ringNode2].netConnections.size() <= 3) {
        logger->debug("No valid move: Switch would result in a ring size less than 3");
        return false;
    }

    // If the ringNodes are a member of fixedRings, the move will be allowed if the resulting ring size
    // is within +/- 1 of the original ring size. This is to allow 3/4 membered rings adjacent to the fixedRing
    // to be able to escape being so. This would otherwise be impossible to remove 3/4 membered rings adjacent to a fixedRing.
    if (fixedRings.count(ringNode1) > 0) {
        int currentSize = networkB.nodes[ringNode1].netConnections.size();
        if (currentSize == fixedRings[ringNode1] - 1) {
            logger->debug("No valid move: Switch would violate fixed ring size");
            return false;
        }
    }

    if (fixedRings.count(ringNode2) > 0) {
        int currentSize = networkB.nodes[ringNode2].netConnections.size();
        if (currentSize == fixedRings[ringNode2] - 1) {
            logger->debug("No valid move: Switch would violate fixed ring size");
            return false;
        }
    }

    if (fixedRings.count(ringNode3) > 0) {
        int currentSize = networkB.nodes[ringNode3].netConnections.size();
        if (currentSize == fixedRings[ringNode3] + 1) {
            logger->debug("No valid move: Switch would violate fixed ring size");
            return false;
        }
    }

    if (fixedRings.count(ringNode4) > 0) {
        int currentSize = networkB.nodes[ringNode4].netConnections.size();
        if (currentSize == fixedRings[ringNode4] + 1) {
            logger->debug("No valid move: Switch would violate fixed ring size");
            return false;
        }
    }

    // check move will not violate dual connectivity limits
    if (networkB.nodes[ringNode1].netConnections.size() == minRingSize ||
        networkB.nodes[ringNode2].netConnections.size() == minRingSize ||
        networkB.nodes[ringNode3].netConnections.size() == maxRingSize ||
        networkB.nodes[ringNode4].netConnections.size() == maxRingSize) {
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
        logger->debug("            |       |                                \\ {:03}  / ", ringNode2);
        logger->debug("            |  {:03}  |                             \\    /", ringNode2);
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
        logger->debug("            |  {:03}  |                                /   \\, ringNode1");
        logger->debug("            |       |                                 / {:03} \\ ", ringNode1);
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
        logger->debug("               {:03}                                      {:03}", baseNode9, baseNode9);
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

/**
 * @brief Find a common base node connection between a base node and ring node excluding a given node
 * @param baseNode ID of the base node
 * @param ringNode ID of the ring node
 * @param excludeNode ID of the node to be excluded
 * @return ID of the common connection
 * @throw std::runtime_error if the associated node cannot be found
 */
int LinkedNetwork::findCommonConnection(const int &baseNode, const int &ringNode, const int &excludeNode) const {
    // Find node that shares baseNode and ringNode but is not excludeNode
    std::unordered_set<int> commonConnections = intersectVectors(networkA.nodes[baseNode].netConnections, networkB.nodes[ringNode].dualConnections);
    commonConnections.erase(excludeNode);
    if (commonConnections.size() != 1) {
        throw std::runtime_error("Could not find common base node for base node " + std::to_string(baseNode) +
                                 " and ring node " + std::to_string(ringNode) +
                                 " excluding node " + std::to_string(excludeNode));
    }
    return *commonConnections.begin();
}
/**
 * @brief Find a common ring connection between two base nodes that exlcudes a given node
 * @param baseNode1 ID of the first base node
 * @param baseNode2 ID of the second base node
 * @param excludeNode ID of the ring node to be excluded
 * @return ID of the common ring connection
 * @throw std::runtime_error if the associated node cannot be found
 */
int LinkedNetwork::findCommonRing(const int &baseNode1, const int &baseNode2, const int &excludeNode) const {
    // Find node that shares baseNode1 and baseNode2 but is not excludeNode
    std::unordered_set<int> commonRings = intersectVectors(networkA.nodes[baseNode1].dualConnections, networkA.nodes[baseNode2].dualConnections);
    commonRings.erase(excludeNode);
    if (commonRings.size() != 1) {
        throw std::runtime_error("Could not find common ring node for base node " + std::to_string(baseNode1) +
                                 " and base node " + std::to_string(baseNode2) +
                                 " excluding ring node " + std::to_string(excludeNode));
    }
    return *commonRings.begin();
}

/**
 * @brief Switch the BSS network by breaking and making bonds
 * @param bondBreaks the bonds to break (vector of pairs)
 * @param bondMakes the bonds to make (vector of pairs)
 */
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

    auto &nodeA1 = networkA.nodes[atom1];
    auto &nodeA2 = networkA.nodes[atom2];
    auto &nodeA4 = networkA.nodes[atom4];
    auto &nodeA5 = networkA.nodes[atom5];

    auto &nodeB1 = networkB.nodes[ringNode1];
    auto &nodeB2 = networkB.nodes[ringNode2];
    auto &nodeB3 = networkB.nodes[ringNode3];
    auto &nodeB4 = networkB.nodes[ringNode4];

    // A-A connectivities
    replaceValue(nodeA1.netConnections, atom5, atom4);
    replaceValue(nodeA2.netConnections, atom4, atom5);
    replaceValue(nodeA4.netConnections, atom2, atom1);
    replaceValue(nodeA5.netConnections, atom1, atom2);

    // A-B connectvities
    replaceValue(nodeA1.dualConnections, ringNode1, ringNode4);
    replaceValue(nodeA2.dualConnections, ringNode2, ringNode3);

    // B-B connectivities
    deleteByValues(nodeB1.netConnections, ringNode2);
    deleteByValues(nodeB2.netConnections, ringNode1);
    nodeB3.netConnections.emplace_back(ringNode4);
    nodeB4.netConnections.emplace_back(ringNode3);

    // B-A connectivities
    deleteByValues(nodeB1.dualConnections, atom1);
    deleteByValues(nodeB2.dualConnections, atom2);

    nodeB3.dualConnections.emplace_back(atom2);
    nodeB4.dualConnections.emplace_back(atom1);
}

/**
 * @brief Restores the initial state of the network by assigning Node objects to their initial states
 * @param initialInvolvedNodesA the initial state of the nodes in network A
 * @param initialInvolvedNodesB the initial state of the nodes in network B
 */
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
 * @brief Sets the coordinates of nodes in the network. Ring network coordinates are recalculated.
 * @param coords Coordinates of the base nodes
 * @throw std::invalid_argument if the number of coordinates does not match the number of nodes in network A
 */
void LinkedNetwork::pushCoords(const std::vector<double> &coords) {
    if (coords.size() != 2 * networkA.nodes.size()) {
        throw std::invalid_argument("Number of coordinates does not match number of nodes in network A");
    }

    // Sync A coordinates
    for (int i = 0; i < networkA.nodes.size(); ++i) {
        networkA.nodes[i].crd[0] = coords[2 * i];
        networkA.nodes[i].crd[1] = coords[2 * i + 1];
    }
    // Centre all the rings relative to their dual connections
    networkB.centreRings(networkA);
}

/**
 * @brief Checks if network is consistent and logs any inconsistencies
 * @return true if all connectivities are reciprocated and all base nodes have clockwise neighbours
 * and ring sizes that are not neighbours of fixedRings are within range, false otherwise
 */
bool LinkedNetwork::checkConsistency() {
    logger->info("Checking consistency...");
    bool consistent = true;
    std::for_each(networkA.nodes.begin(), networkA.nodes.end(), [this, &consistent](const Node &node) {
        std::for_each(node.netConnections.begin(), node.netConnections.end(), [this, &node, &consistent](const int &cnx) {
            if (!vectorContains(networkA.nodes[cnx].netConnections, node.id)) {
                logger->error("Node {} base has neighbour {} but neighbour does not have node as neighbour", node.id, cnx);
                consistent = false;
            }
        });
    });
    std::unordered_set<int> fixedRingNeighbours = {};
    std::for_each(fixedRings.begin(), fixedRings.end(), [this, &fixedRingNeighbours](const std::pair<int, int> &ring) {
        fixedRingNeighbours.insert(networkB.nodes[ring.first].netConnections.begin(), networkB.nodes[ring.first].netConnections.end());
    });
    std::for_each(networkB.nodes.begin(), networkB.nodes.end(), [this, &consistent, &fixedRingNeighbours](const Node &node) {
        std::for_each(node.netConnections.begin(), node.netConnections.end(), [this, &node, &consistent](const int &cnx) {
            if (!vectorContains(networkB.nodes[cnx].netConnections, node.id)) {
                logger->error("Node {} ring has neighbour {} but neighbour does not have node as neighbour", node.id, cnx);
                consistent = false;
            }
        });
        if (fixedRingNeighbours.count(node.id) == 0 && (node.netConnections.size() > maxRingSize || node.netConnections.size() < minRingSize)) {
            logger->error("Node {} ring has {} neighbours, which is outside the allowed range", node.id, node.netConnections.size());
            consistent = false;
        }
    });
    std::for_each(networkA.nodes.begin(), networkA.nodes.end(), [this, &consistent](Node &node) {
        std::for_each(node.dualConnections.begin(), node.dualConnections.end(), [this, &node, &consistent](const int &cnx) {
            if (!vectorContains(networkB.nodes[cnx].dualConnections, node.id)) {
                logger->error("Node {} base has ring neighbour {} but ring neighbour does not have node as ring neighbour", node.id, cnx);
                consistent = false;
            }
        });
    });
    std::for_each(networkB.nodes.begin(), networkB.nodes.end(), [this, &consistent](const Node &node) {
        std::for_each(node.dualConnections.begin(), node.dualConnections.end(), [this, &node, &consistent](const int &cnx) {
            if (!vectorContains(networkA.nodes[cnx].dualConnections, node.id)) {
                logger->error("Node {} ring has ring neighbour {} but ring neighbour does not have node as ring neighbour", node.id, cnx);
                consistent = false;
            }
        });
    });
    return checkAllClockwiseNeighbours() && consistent;
}

/**
 * @brief Wraps coordinates out of bounds back into the periodic box, only for 2 dimensions
 * @param coords Coordinates to be wrapped (1D vector of pairs)
 */
void LinkedNetwork::wrapCoords(std::vector<double> &coords) const {
    int dim = 0;
    for (auto &coord : coords) {
        coord = fmod(coord, dimensions[dim]);
        dim = 1 - dim; // Alternate between 0 and 1
    }
}

/**
 * @brief Writes the network to files
*/
void LinkedNetwork::write() const {
    networkA.write();
    networkB.write();
}

/**
 * @brief Checks if a given node has clockwise neighbours
 * @param nodeID ID of the node to check
 * @return true if the node has clockwise neighbours, false otherwise
 */
bool LinkedNetwork::checkClockwiseNeighbours(const int &nodeID) const {
    Node node = networkA.nodes[nodeID];
    double prevAngle = getClockwiseAngle(node.crd, networkA.nodes[node.netConnections.back()].crd, dimensions);

    // The angle a neighbour makes with the x axis should always increase, EXCEPT when the angle wraps around from 2π to 0
    // This should only happen once if the neighbours are clockwise
    int timesDecreased = 0;
    for (int id : node.netConnections) {
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
    const int lastNeighbourCoordsID = networkA.nodes[nodeID].netConnections.back();
    const std::vector<double> lastNeighbourCoord = {coords[2 * lastNeighbourCoordsID], coords[2 * lastNeighbourCoordsID + 1]};
    double prevAngle = getClockwiseAngle(nodeCoord, lastNeighbourCoord, dimensions);
    int timesDecreased = 0;
    for (const auto &id : networkA.nodes[nodeID].netConnections) {
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
 * @return true if all nodes have clockwise neighbours, false otherwise
 */
bool LinkedNetwork::checkAllClockwiseNeighbours() const {
    bool allClockwise = true;
    for (int nodeID = 0; nodeID < networkA.nodes.size(); ++nodeID) {
        if (!checkClockwiseNeighbours(nodeID)) {
            logger->warn("Node {} has anticlockwise neighbours: {}", nodeID, vectorToString(networkA.nodes[nodeID].netConnections));
            allClockwise = false;
        }
    }
    return allClockwise;
}

void LinkedNetwork::arrangeNeighboursClockwise(const int &nodeID, const std::vector<double> &coords) {
    // Get the coordinates of the center node
    const std::vector<double> nodeCoord = {coords[2 * nodeID], coords[2 * nodeID + 1]};

    // Get the neighbour IDs of the center node
    std::vector<int> &neighbours = networkA.nodes[nodeID].netConnections;

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
    for (size_t i = 0; i < neighbours.size(); ++i) {
        neighbours[i] = neighbourAngles[i].first;
    }
}

void LinkedNetwork::arrangeNeighboursClockwise(const std::unordered_set<int> &nodeIDs, const std::vector<double> &coords) {
    std::for_each(nodeIDs.begin(), nodeIDs.end(), [this, &coords](int nodeID) {
        arrangeNeighboursClockwise(nodeID, coords);
    });
}

/**
 * @brief Checks if all angles around all nodes are less than or equal to maximum angle
 * @param coords Coordinates of all nodes as a 1D vector of coordinate pairs
 * @return true if all angles are convex, false otherwise
 */
bool LinkedNetwork::checkAnglesWithinRange(const std::vector<double> &coords) {
    for (int nodeID = 0; nodeID < networkA.nodes.size(); ++nodeID) {
        arrangeNeighboursClockwise(nodeID, coords);
        for (int i = 0; i < networkA.nodes[nodeID].netConnections.size(); ++i) {
            int neighbourID = networkA.nodes[nodeID].netConnections[i];
            int nextNeighbourID = networkA.nodes[nodeID].netConnections[(i + 1) % networkA.nodes[nodeID].netConnections.size()];
            std::vector<double> v1 = pbcVector(std::vector<double>{coords[neighbourID * 2], coords[neighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            std::vector<double> v2 = pbcVector(std::vector<double>{coords[nextNeighbourID * 2], coords[nextNeighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            double angle = getClockwiseAngleBetweenVectors(v1, v2);
            if (angle > maximumAngle) {
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
 * @return true if all angles are within range, false otherwise
 */
bool LinkedNetwork::checkAnglesWithinRange(const std::unordered_set<int> &nodeIDs, const std::vector<double> &coords) {
    for (int nodeID : nodeIDs) {
        arrangeNeighboursClockwise(nodeID, coords);
        for (int i = 0; i < networkA.nodes[nodeID].netConnections.size(); ++i) {
            int neighbourID = networkA.nodes[nodeID].netConnections[i];
            int nextNeighbourID = networkA.nodes[nodeID].netConnections[(i + 1) % networkA.nodes[nodeID].netConnections.size()];
            std::vector<double> v1 = pbcVector(std::vector<double>{coords[neighbourID * 2], coords[neighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            std::vector<double> v2 = pbcVector(std::vector<double>{coords[nextNeighbourID * 2], coords[nextNeighbourID * 2 + 1]},
                                               std::vector<double>{coords[nodeID * 2], coords[nodeID * 2 + 1]}, dimensions);
            double angle = getClockwiseAngleBetweenVectors(v1, v2);
            if (angle > maximumAngle) {
                return false;
            }
        }
    }
    return true;
}

/**
 * @brief Checks if the bonds to a given node are within the maximum bond length
 * @param nodeID ID of the node to check
 * @param coords Coordinates of all nodes as a 1D vector of coordinate pairs
 * @return true if all bonds are within the maximum bond length, false otherwise
 */
bool LinkedNetwork::checkBondLengths(const int &nodeID, const std::vector<double> &coords) const {
    for (const auto &neighbourID : networkA.nodes[nodeID].netConnections) {
        std::vector<double> pbcVec = pbcVector(std::vector<double>{coords[2 * nodeID], coords[2 * nodeID + 1]},
                                               std::vector<double>{coords[2 * neighbourID], coords[2 * neighbourID + 1]}, dimensions);
        if (std::sqrt(pbcVec[0] * pbcVec[0] + pbcVec[1] * pbcVec[1]) > maximumBondLength) {
            return false;
        }
    }
    return true;
}

/**
 * @brief Checks if the bonds to the given nodes are within the maximum bond length
 * @param nodeIDs IDs of the nodes to check
 * @param coords Coordinates of all nodes as a 1D vector of coordinate pairs
 * @return true if all bonds are within the maximum bond length, false otherwise
 */
bool LinkedNetwork::checkBondLengths(const std::unordered_set<int> &nodeIDs, const std::vector<double> &coords) const {
    return std::all_of(nodeIDs.begin(), nodeIDs.end(), [this, &coords](int nodeID) {
        return checkBondLengths(nodeID, coords);
    });
}

/**
 * @brief Gets the order of a given array of 4 ring IDs using the centre of all of their coordinates
 * @param ringNodeIDs IDs of the rings to get the order of
 * @return clockwise or anti-clockwise
 */
Direction LinkedNetwork::getRingsDirection(const std::vector<int> &ringNodeIDs) const {
    if (ringNodeIDs.size() != 4) {
        throw std::invalid_argument("Error getting ring direction, ringNodeIDs size is not 4");
    }
    std::vector<double> midCoords(2, 0.0);
    for (int id : ringNodeIDs) {
        midCoords[0] += currentCoords[id * 2];
        midCoords[1] += currentCoords[id * 2 + 1];
    }
    midCoords[0] /= 4;
    midCoords[1] /= 4;
    int timesDecreased = 0;
    double prevAngle = getClockwiseAngle(midCoords,
                                         {currentCoords[ringNodeIDs.back() * 2], currentCoords[ringNodeIDs.back() * 2 + 1]},
                                         dimensions);
    for (int id : ringNodeIDs) {
        double angle = getClockwiseAngle(midCoords, {currentCoords[id * 2], currentCoords[id * 2 + 1]}, dimensions);
        if (angle < prevAngle) {
            timesDecreased++;
            if (timesDecreased == 2) {
                logger->debug("Anticlockwise");
                return Direction::ANTICLOCKWISE;
            }
        }
        prevAngle = angle;
    }
    logger->debug("Clockwise");
    return Direction::CLOCKWISE;
}

/**
 * @brief Calculates the new coordinates of a pair of atoms if they were rotated by 90 degrees in the given direction
 * @param atomID1 ID of the first atom in the bond
 * @param atomID2 ID of the second atom in the bond
 * @param direct Direction to rotate the bond
 * @return Pair of vectors containing the new coordinates of the atoms
 */
std::tuple<std::vector<double>, std::vector<double>> LinkedNetwork::rotateBond(const int &atomID1, const int &atomID2,
                                                                               const Direction &direct) const {
    logger->debug("Rotating bond between atoms {} and {}", atomID1, atomID2);
    std::vector<double> atom1Coord = {currentCoords[atomID1 * 2], currentCoords[atomID1 * 2 + 1]};
    std::vector<double> atom2Coord = {currentCoords[atomID2 * 2], currentCoords[atomID2 * 2 + 1]};

    // Calculate the center point
    double centerX = (atom1Coord[0] + atom2Coord[0]) / 2.0;
    double centerY = (atom1Coord[1] + atom2Coord[1]) / 2.0;

    // Translate the points to the origin
    atom1Coord[0] -= centerX;
    atom1Coord[1] -= centerY;
    atom2Coord[0] -= centerX;
    atom2Coord[1] -= centerY;

    std::vector<double> rotatedAtom1Coord;
    std::vector<double> rotatedAtom2Coord;
    // Rotate the points by 90 degrees
    if (direct == Direction::CLOCKWISE) {
        rotatedAtom1Coord = {atom1Coord[1], -atom1Coord[0]};
        rotatedAtom2Coord = {atom2Coord[1], -atom2Coord[0]};
    } else {
        rotatedAtom1Coord = {-atom1Coord[1], atom1Coord[0]};
        rotatedAtom2Coord = {-atom2Coord[1], atom2Coord[0]};
    }
    // Translate the points back
    rotatedAtom1Coord[0] += centerX;
    rotatedAtom1Coord[1] += centerY;
    rotatedAtom2Coord[0] += centerX;
    rotatedAtom2Coord[1] += centerY;

    wrapCoords(rotatedAtom1Coord);
    wrapCoords(rotatedAtom2Coord);

    return {rotatedAtom1Coord, rotatedAtom2Coord};
}

/**
 * @brief Gets the areas of all rings in the network
 * @return Vector of areas of all rings
 */
std::vector<double> LinkedNetwork::getRingAreas() const {
    std::vector<double> ringAreas;
    ringAreas.reserve(networkB.nodes.size());
    for (size_t ringNodeID = 0; ringNodeID < networkB.nodes.size(); ringNodeID++) {
        std::vector<std::vector<double>> baseNodeCoords;
        baseNodeCoords.reserve(networkB.nodes[ringNodeID].dualConnections.size());
        for (size_t baseNodeID : networkB.nodes[ringNodeID].dualConnections) {
            std::vector<double> pbcVec = pbcVector({currentCoords[baseNodeID * 2], currentCoords[baseNodeID * 2 + 1]},
                                                   {currentCoords[ringNodeID * 2], currentCoords[ringNodeID * 2 + 1]}, dimensions);
            baseNodeCoords.push_back(pbcVec);
        }
        double ringArea = calculatePolygonArea(baseNodeCoords);
        ringAreas.push_back(ringArea);
    }
    return ringAreas;
}