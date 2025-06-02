#include "linked_network.h"
#include "file_tools.h"
#include "random_number_generator.h"
#include "switch_move.h"
#include "types.h"
#include "vector_tools.h"
#include <algorithm>
#include <array>
#include <chrono>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <format>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <ranges>
#include <spdlog/spdlog.h>
#include <sstream>

/**
 * @brief Construct by loading networks from files
 * @param inputData the input data object
 * @param loggerArg the logger object
 */
LinkedNetwork::LinkedNetwork(const InputData &inputData,
                             const LoggerPtr &loggerArg)
    : minRingSize(inputData.minRingSize), maxRingSize(inputData.maxRingSize),
      selectionType(inputData.randomOrWeighted),
      metropolisCondition(Metropolis()), weightedDecay(inputData.weightedDecay),
      maximumBondLength(inputData.maximumBondLength),
      maximumAngle(inputData.maximumAngle * M_PI / 180),
      writeMovie(inputData.writeMovie), logger(loggerArg) {
  networkA = Network(NetworkType::BASE_NETWORK, logger);
  networkB = Network(NetworkType::DUAL_NETWORK, logger);
  if (inputData.isFixRingsEnabled) {
    findFixedRings(std::filesystem::path("./input_files") / "bss_network" /
                   "fixed_rings.txt");
    findFixedNodes();
  } else {
    logger->info("Fixed rings disabled, setting number of fixed rings to 0.");
  }
  if (size_t loadedMinRingSize = networkB.getMinConnections();
      loadedMinRingSize < minRingSize) {
    logger->warn("Loaded network has a min ring size of {} which is lower than "
                 "input file's {}",
                 loadedMinRingSize, minRingSize);
  }
  if (size_t loadedMaxRingSize = networkB.getMaxConnections();
      loadedMaxRingSize > maxRingSize) {
    logger->warn("Loaded network has a max ring size of {} which is higher "
                 "than input file's {}",
                 loadedMaxRingSize, maxRingSize);
  }
  dimensions = networkA.dimensions;
  centreCoords = {dimensions[0] / 2, dimensions[1] / 2};

  lammpsManager = LAMMPSManager(logger);
  if (writeMovie) {
    lammpsManager.startMovie();
    lammpsManager.writeMovie();
  }
  lammpsManager.minimiseNetwork();
  energy = lammpsManager.getPotentialEnergy();
  pushCoords(lammpsManager.getCoords());
  weights.resize(networkA.nodes.size());
  updateWeights();
}

/**
 * @brief Read the fixed_rings.txt file and populate fixedRings with integers
 * from each line
 * @param isFixedRingsEnabled boolean to enable or disable fixed rings
 * @param filename the name of the input file
 */
void LinkedNetwork::findFixedRings(const std::string &filePath) {
  const std::unordered_set<uint16_t> fixedRingIDs =
      readFixedRings(filePath, logger);
  std::ranges::for_each(fixedRingIDs, [this](const uint16_t fixedRingID) {
    this->fixedRings[fixedRingID] =
        networkB.nodes[fixedRingID].numConnections();
  });
  logger->info("Number of fixed rings: {}", fixedRings.size());
  logger->info("Fixed ring info (ID: Size): {}", mapToString(fixedRings));
}

/**
 * @brief Populate fixedNodes with all base nodes that are a member of all fixed
 * rings
 */
void LinkedNetwork::findFixedNodes() {
  std::ranges::for_each(
      fixedRings,
      [this](const std::pair<const uint16_t, size_t> &fixedRingPair) {
        this->fixedNodes.insert(
            this->networkB.nodes[fixedRingPair.first].dualConnections.begin(),
            this->networkB.nodes[fixedRingPair.first].dualConnections.end());
      });
}

std::optional<SwitchMove> LinkedNetwork::getSwitchMove() {
  // Number of bonds in a network of N nodes with degree 3 is 3 * N / 2
  const size_t maxAttempts = networkA.nodes.size() * 3 / 2;
  size_t attempts = 0;
  std::set<std::set<uint16_t>> failedConnections;
  while (true) {
    if (attempts >= maxAttempts) {
      logger->error("Cannot find any valid switch moves");
      break;
    }
    std::array<std::array<uint16_t, 2>, 2> result = pickRandomConnection();
    auto &[baseConnection, dualConnection] = result;
    auto &[baseNode1, baseNode2] = baseConnection;
    auto &[ringNode1, ringNode2] = dualConnection;
    if (failedConnections.contains(
            std::set{baseConnection[0], baseConnection[1]})) {
      // We have already tried this connection
      continue;
    }
    logger->debug("Random connection chosen. Base: {}-{} Ring: {}-{}",
                  baseNode1, baseNode2, ringNode1, ringNode2);
    try {
      return std::optional<SwitchMove>{
          this->genSwitchMove(baseNode1, baseNode2, ringNode1, ringNode2)};
    } catch (const SwitchMoveException &e) {
      logger->debug("Invalid move: {}", e.what());
      failedConnections.insert(std::set{baseConnection[0], baseConnection[1]});
      attempts++;
    }
  }
  return std::nullopt;
}

/**
 * @brief Perform a bond switch, evaluate energy, accept or
 * reject
 * @param temperature The temperature of the system
 */
void LinkedNetwork::performBondSwitch(const double temperature) {
  logger->debug("Finding move...");
  const std::optional<SwitchMove> switchMoveOpt = getSwitchMove();
  if (!switchMoveOpt.has_value()) {
    logger->error("Cannot find any valid switch moves");
    throw std::runtime_error("Cannot find any valid switch moves");
  }
  const SwitchMove &switchMove = switchMoveOpt.value();

  numSwitches++;
  logger->debug("Switch number: {}", numSwitches);

  // Save current state
  double initialEnergy = energy;

  std::vector<Node> initialInvolvedNodesA;
  for (const uint16_t &id : switchMove.involvedBaseNodes) {
    initialInvolvedNodesA.push_back(networkA.nodes[id]);
  }
  std::vector<Node> initialInvolvedNodesB;
  for (const std::array<uint16_t, 2> &bond : switchMove.ringBondBreakMake) {
    for (const uint16_t &id : bond) {
      initialInvolvedNodesB.push_back(networkB.nodes[id]);
    }
  }

  // Switch and geometry optimise
  logger->debug("Switching BSS Network...");
  applyMove(switchMove.bondBreaks, switchMove.ringBondBreakMake);
  std::array<double, 2> rotatedCoord1;
  std::array<double, 2> rotatedCoord2;
  std::array<uint16_t, 4> orderedRingNodes = {
      switchMove.ringBondBreakMake[0][1], switchMove.ringBondBreakMake[1][1],
      switchMove.ringBondBreakMake[0][0], switchMove.ringBondBreakMake[1][0]};
  std::tie(rotatedCoord1, rotatedCoord2) = this->networkA.getRotatedBond(
      switchMove.selectedBaseBond, getRingsDirection(orderedRingNodes));
  logger->debug("Switching LAMMPS Network...");
  lammpsManager.performSwitch(switchMove, rotatedCoord1, rotatedCoord2);

  logger->debug("Minimising network...");
  lammpsManager.minimiseNetwork();

  logger->debug("Accepting or rejecting...");
  std::vector<std::array<double, 2>> potentialCoords =
      lammpsManager.getCoords();
  // TODO - Check angles are within range
  if (!checkBondLengths(switchMove.involvedBaseNodes, potentialCoords)) {
    logger->debug("Rejected move: bond lengths are not within range");
    failedBondLengthChecks++;
    rejectMove(initialInvolvedNodesA, initialInvolvedNodesB, switchMove);
    return;
  }
  if (!this->checkAngles(switchMove.involvedBaseNodes, potentialCoords)) {
    logger->debug("Rejected move: angles are not within range");
    failedAngleChecks++;
    rejectMove(initialInvolvedNodesA, initialInvolvedNodesB, switchMove);
    return;
  }
  double finalEnergy = lammpsManager.getPotentialEnergy();
  if (std::isnan(finalEnergy)) {
    throw std::runtime_error("Final energy is NaN");
  }
  if (!metropolisCondition.acceptanceCriterion(finalEnergy, initialEnergy,
                                               temperature)) {
    logger->debug("Rejected move: failed Metropolis criterion: Ei = {:.3f} Eh, "
                  "Ef = {:.3f} Eh",
                  initialEnergy, finalEnergy);
    failedEnergyChecks++;
    rejectMove(initialInvolvedNodesA, initialInvolvedNodesB, switchMove);
    return;
  }

  this->acceptMove(potentialCoords, initialEnergy, finalEnergy);
}

void LinkedNetwork::acceptMove(
    const std::vector<std::array<double, 2>> &newCoords,
    const double initialEnergy, const double finalEnergy) {
  logger->debug("Accepted Move: Ei = {:.3f} Eh, Ef = {:.3f} Eh", initialEnergy,
                finalEnergy);
  this->numAcceptedSwitches++;
  logger->debug("Syncing LAMMPS coordinates to BSS coordinates...");
  this->pushCoords(newCoords);
  this->updateWeights();
  this->energy = finalEnergy;
  if (this->writeMovie) {
    this->lammpsManager.writeMovie();
  }
}

void LinkedNetwork::rejectMove(const std::vector<Node> &initialInvolvedNodesA,
                               const std::vector<Node> &initialInvolvedNodesB,
                               const SwitchMove &switchMove) {
  logger->debug("Reverting BSS Network...");
  this->revertMove(initialInvolvedNodesA, initialInvolvedNodesB);
  logger->debug("Reverting LAMMPS Network...");
  this->lammpsManager.revertSwitch(switchMove);
  logger->debug("Syncing LAMMPS coordinates to BSS coordinates...");
  this->lammpsManager.setCoords(networkA.getCoords());
}

void LinkedNetwork::showCoords(
    const std::vector<std::array<double, 2>> &coords) const {
  std::ranges::for_each(coords, [this](const std::array<double, 2> &coord) {
    logger->debug("{} {}", coord[0], coord[1]);
  });
}

/**
 * @brief Scale the entire linked network by a given factor
 * @param scaleFactor the factor to scale the network by
 */
void LinkedNetwork::rescale(double scaleFactor) {
  containerMultiply(dimensions, scaleFactor);
  networkA.rescale(scaleFactor);
  networkB.rescale(scaleFactor);
}

void LinkedNetwork::updateWeights() {
  if (selectionType == SelectionType::EXPONENTIAL_DECAY) {
    double boxLength = dimensions[0];
    for (int i = 0; i < networkA.nodes.size(); ++i) {
      double distance =
          networkA.nodes[i].distanceFrom(centreCoords) / boxLength;
      weights[i] = std::exp(-distance * weightedDecay);
    }

    // Normalize weights
    double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    for (double &weight : weights) {
      weight /= total;
    }
  } else { // SelectionType::RANDOM
    std::ranges::fill(weights, 1.0);
  }
}

/**
 * @brief Chooses a random bond in the network and returns the IDs in the bond,
 * the two rings either side and the connection type
 * @return A tuple containing the IDs of the two nodes in the bond, the two
 * rings either side and the connection type
 * @throw std::runtime_error if the nodes in the random connection have
 * coordinations other than 3 or 4.
 */
std::array<std::array<uint16_t, 2>, 2> LinkedNetwork::pickRandomConnection() {
  const Node &randomNode = this->networkA.getRandomNode();
  const Node &randomNodeConnection =
      this->networkA.getRandomNodeConnection(randomNode);

  std::unordered_set<uint16_t> commonRingIDs = intersectSets(
      randomNode.dualConnections, randomNodeConnection.dualConnections);
  // Two connected nodes should always share two ring nodes.
  if (commonRingIDs.size() != 2) {
    throw std::runtime_error(std::format(
        "Selected random connection does not share two ring nodes: {} {}",
        randomNode.id, randomNodeConnection.id));
  }
  // Randomly assign ringNode1 and ringNode2 to those two common ring
  // nodes
  uint16_t sharedRingNode1 =
      *RandomNumberGenerator::getInstance().getRandomElement(
          commonRingIDs.begin(), commonRingIDs.end());
  commonRingIDs.erase(sharedRingNode1);
  uint16_t sharedRingNode2 = *commonRingIDs.begin();

  return {{{randomNode.id, randomNodeConnection.id},
           {sharedRingNode1, sharedRingNode2}}};
}

SwitchMove LinkedNetwork::genSwitchMove(const uint16_t baseNode1,
                                        const uint16_t baseNode2,
                                        const uint16_t ringNode1,
                                        const uint16_t ringNode2) const {
  if (baseNode1 == baseNode2 || ringNode1 == ringNode2) {
    throw SwitchMoveException(std::format(
        "Cannot switch common base nodes: {}-{} or ring nodes: {}-{}",
        baseNode1, baseNode2, ringNode1, ringNode2));
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
  uint16_t baseNode5 =
      this->findCommonConnection(baseNode1, ringNode1, baseNode2);
  uint16_t baseNode6 =
      this->findCommonConnection(baseNode2, ringNode1, baseNode1);
  uint16_t baseNode3 =
      this->findCommonConnection(baseNode1, ringNode2, baseNode2);
  uint16_t baseNode4 =
      this->findCommonConnection(baseNode2, ringNode2, baseNode1);

  uint16_t ringNode3 = this->findCommonRing(baseNode1, baseNode5, ringNode1);
  uint16_t ringNode4 = this->findCommonRing(baseNode2, baseNode6, ringNode1);
  uint16_t baseNode11 =
      this->findCommonConnection(baseNode3, ringNode3, baseNode1);
  uint16_t baseNode7 =
      this->findCommonConnection(baseNode3, ringNode2, baseNode1);
  uint16_t baseNode8 =
      this->findCommonConnection(baseNode4, ringNode2, baseNode2);
  uint16_t baseNode12 =
      this->findCommonConnection(baseNode4, ringNode4, baseNode2);
  uint16_t baseNode14 =
      this->findCommonConnection(baseNode6, ringNode4, baseNode2);
  uint16_t baseNode10 =
      this->findCommonConnection(baseNode6, ringNode1, baseNode2);
  uint16_t baseNode9 =
      this->findCommonConnection(baseNode5, ringNode1, baseNode1);
  uint16_t baseNode13 =
      this->findCommonConnection(baseNode5, ringNode3, baseNode1);

  if (baseNode5 == baseNode6 || baseNode3 == baseNode4) {
    throw SwitchMoveException(
        "No valid move: Selected nodes describe an edge of two edge "
        "sharing triangles");
  }

  if (this->networkB.nodes[ringNode1].numConnections() <= 3 ||
      this->networkB.nodes[ringNode2].numConnections() <= 3) {
    throw SwitchMoveException(
        "No valid move: Switch would result in a ring size less than 3");
  }

  if (!this->checkConnectivityLimits(ringNode1, ringNode2, ringNode3,
                                     ringNode4)) {
    throw SwitchMoveException(
        "No valid move: Switch would violate ring size limits");
  }

  const std::array<std::array<uint16_t, 2>, 2> bondBreaks = {
      {{baseNode1, baseNode5}, {baseNode2, baseNode4}}};
  const std::array<std::array<uint16_t, 2>, 2> bondMakes = {
      {{baseNode1, baseNode4}, {baseNode2, baseNode5}}};
  //                          Break                   Make
  const std::array<std::array<uint16_t, 2>, 2> ringBondBreakMake = {
      {{ringNode1, ringNode2}, {ringNode3, ringNode4}}};

  const std::unordered_set<uint16_t> involvedBaseNodes = {
      baseNode1,  baseNode2,  baseNode3,  baseNode4, baseNode5,
      baseNode6,  baseNode7,  baseNode8,  baseNode9, baseNode10,
      baseNode11, baseNode12, baseNode13, baseNode14};

  std::array<std::array<uint16_t, 3>, 8> angleMakes;
  std::array<std::array<uint16_t, 3>, 8> angleBreaks;
  if (baseNode7 == baseNode4) {
    // 4 membered ringNode2
    angleBreaks[0] = {baseNode3, baseNode4, baseNode2};
    angleBreaks[1] = {baseNode12, baseNode4, baseNode2};
    angleBreaks[2] = {baseNode1, baseNode2, baseNode4};
    angleBreaks[3] = {baseNode6, baseNode2, baseNode4};

    angleMakes[0] = {baseNode3, baseNode4, baseNode1};
    angleMakes[1] = {baseNode12, baseNode4, baseNode1};
    angleMakes[2] = {baseNode3, baseNode1, baseNode4};
    angleMakes[3] = {baseNode2, baseNode1, baseNode4};

    logger->debug(" {:03}------{:03}------{:03}------{:03}             "
                  "{:03}-----{:03}-----{:03}-----{:03} ",
                  baseNode11, baseNode3, baseNode4, baseNode12, baseNode11,
                  baseNode3, baseNode4, baseNode12);
    logger->debug(
        "            |       |                                \\ {:03}  / ",
        ringNode2);
    logger->debug("            |  {:03}  |                             \\    /",
                  ringNode2);
    logger->debug(
        "            |       |                                  {:03}",
        baseNode1);
  } else if (baseNode7 == baseNode8) {
    // 5 membered ringNode2
    angleBreaks = {baseNode7, baseNode4, baseNode2, baseNode12,
                   baseNode4, baseNode2, baseNode1, baseNode2,
                   baseNode4, baseNode6, baseNode2, baseNode4};
    angleMakes = {baseNode7, baseNode4, baseNode1, baseNode12,
                  baseNode4, baseNode1, baseNode2, baseNode1,
                  baseNode4, baseNode3, baseNode1, baseNode4};
    logger->debug("");
    logger->debug(
        "               {:03}                                      {:03}",
        baseNode7, baseNode7);
    logger->debug(
        "            /      \\                                   /   \\");
    logger->debug(
        "           /        \\                                 /     \\");
    logger->debug(" {:03}-----{:03}   {:03}  {:03}-----{:03}             "
                  "{:03}-----{:03} {:03} {:03}-----{:03}",
                  baseNode11, baseNode3, ringNode2, baseNode4, baseNode12,
                  baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
    logger->debug(
        "          \\          /                                \\     /");
    logger->debug(
        "           \\        /                                   {:03}",
        baseNode1);
  } else {
    // 6+ membered ringNode2
    angleBreaks = {baseNode2, baseNode4,  baseNode8, baseNode2,
                   baseNode4, baseNode12, baseNode4, baseNode2,
                   baseNode1, baseNode4,  baseNode2, baseNode6};
    angleMakes = {baseNode1, baseNode4,  baseNode8, baseNode1,
                  baseNode4, baseNode12, baseNode4, baseNode1,
                  baseNode2, baseNode4,  baseNode1, baseNode3};
    logger->debug("");
    logger->debug("           {:03}~~~~~{:03}                              "
                  "{:03}~~~~~{:03}",
                  baseNode7, baseNode8, baseNode7, baseNode8);
    logger->debug(
        "           /        \\                                |       |");
    logger->debug("          /          \\                      "
                  "{:03}-----{:03} {:03} {:03}-----{:03}",
                  baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
    logger->debug(" {:03}-----{:03}   {:03}   {:03}-----{:03}                  "
                  "    \\     /",
                  baseNode11, baseNode3, ringNode2, baseNode4, baseNode12);
    logger->debug(
        "          \\          /                                 \\   /");
    logger->debug(
        "           \\        /                                   {:03}",
        baseNode1);
  }
  logger->debug("    {:03}    {:03}-----{:03}   {:03}          --->     {:03}  "
                "     |      {:03}",
                ringNode3, baseNode1, baseNode2, ringNode4, ringNode3,
                ringNode4);
  if (baseNode5 == baseNode10) {
    // 4 membered ringNode1
    angleMakes[4] = {baseNode13, baseNode5, baseNode2};
    angleMakes[5] = {baseNode6, baseNode5, baseNode2};
    angleMakes[6] = {baseNode1, baseNode2, baseNode5};
    angleMakes[7] = {baseNode6, baseNode2, baseNode5};
    logger->debug(
        "            |       |                                   {:03}",
        baseNode2);
    logger->debug("            |  {:03}  |                                /   "
                  "\\, ringNode1");
    logger->debug(
        "            |       |                                 / {:03} \\ ",
        ringNode1);
    logger->debug(" {:03}-------{:03}-----{:03}-------{:03}            "
                  "{:03}-----{:03}-----{:03}-----{:03} ",
                  baseNode13, baseNode5, baseNode6, baseNode14, baseNode13,
                  baseNode5, baseNode6, baseNode14);
    logger->debug("");
  } else if (baseNode9 == baseNode10) {
    // 5 membered ringNode1
    angleBreaks[4] = {baseNode13, baseNode5, baseNode1};
    angleBreaks[5] = {baseNode9, baseNode5, baseNode1};
    angleBreaks[6] = {baseNode3, baseNode1, baseNode5};
    angleBreaks[7] = {baseNode2, baseNode1, baseNode5};

    angleMakes[4] = {baseNode13, baseNode5, baseNode2};
    angleMakes[5] = {baseNode9, baseNode5, baseNode2};
    angleMakes[6] = {baseNode1, baseNode2, baseNode5};
    angleMakes[7] = {baseNode6, baseNode2, baseNode5};
    logger->debug(
        "           /        \\                                   {:03}",
        baseNode2);
    logger->debug(
        "          /          \\                                /     \\");
    logger->debug(" {:03}-----{:03}   {:03}  {:03}-----{:03}             "
                  "{:03}-----{:03} {:03} {:03}-----{:03}",
                  baseNode13, baseNode5, ringNode1, baseNode6, baseNode14,
                  baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
    logger->debug(
        "           \\        /                                 \\     /");
    logger->debug(
        "            \\      /                                   \\   /");
    logger->debug(
        "               {:03}                                      {:03}",
        baseNode9, baseNode9);
    logger->debug("");
  } else {
    // 6+ membered ringNode1
    angleBreaks[4] = {baseNode1, baseNode5, baseNode9};
    angleBreaks[5] = {baseNode1, baseNode5, baseNode13};
    angleBreaks[6] = {baseNode5, baseNode1, baseNode2};
    angleBreaks[7] = {baseNode5, baseNode1, baseNode3};

    angleMakes[4] = {baseNode2, baseNode5, baseNode9};
    angleMakes[5] = {baseNode2, baseNode5, baseNode13};
    angleMakes[6] = {baseNode1, baseNode2, baseNode5};
    angleMakes[7] = {baseNode6, baseNode2, baseNode5};
    logger->debug(
        "           /        \\                                   {:03}",
        baseNode2);
    logger->debug(
        "          /          \\                                 /   \\");
    logger->debug(" {:03}-----{:03}   {:03}   {:03}-----{:03}                  "
                  "    /     \\",
                  baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
    logger->debug("          \\          /                      "
                  "{:03}-----{:03} {:03} {:03}-----{:03}",
                  baseNode13, baseNode5, ringNode1, baseNode6, baseNode14);
    logger->debug(
        "           \\        /                                |       |");
    logger->debug("           {:03}~~~~~{:03}                              "
                  "{:03}~~~~~{:03}",
                  baseNode9, baseNode10, baseNode9, baseNode10);
    logger->debug("");
  }
  return SwitchMove({{baseNode1, baseNode2}}, {{ringNode1, ringNode2}},
                    bondBreaks, bondMakes, angleBreaks, angleMakes,
                    ringBondBreakMake, involvedBaseNodes);
}

bool LinkedNetwork::checkConnectivityLimits(const uint16_t ringNode1,
                                            const uint16_t ringNode2,
                                            const uint16_t ringNode3,
                                            const uint16_t ringNode4) const {

  if (this->networkB.nodes[ringNode1].numConnections() <= 3 ||
      this->networkB.nodes[ringNode2].numConnections() <= 3) {
    logger->debug(
        "No valid move: Switch would result in a ring size less than 3");
    return false;
  }

  // If the ringNodes are a member of fixedRings, the move will be allowed if
  // the resulting ring size is within +/- 1 of the original ring size. This is
  // to allow 3/4 membered rings adjacent to the fixedRing to be able to escape
  // being so. This would otherwise be impossible to remove 3/4 membered rings
  // adjacent to a fixedRing.
  for (const uint16_t &ringNode : {ringNode1, ringNode2}) {
    if (this->fixedRings.contains(ringNode)) {
      if (size_t currentSize = this->networkB.nodes[ringNode].numConnections();
          currentSize == fixedRings.at(ringNode) - 1) {
        // Ring nodes 1 and 2 will decrease in size, so if statement is -1 here
        logger->debug("No valid move: Switch would violate fixed ring size");
        return false;
      }
    }
  }

  for (const uint16_t &ringNode : {ringNode3, ringNode4}) {
    if (fixedRings.contains(ringNode)) {
      if (size_t currentSize = networkB.nodes[ringNode].numConnections();
          currentSize == fixedRings.at(ringNode) + 1) {
        // Ring nodes 3 and 4 will increase in size, so if statement is +1 here
        logger->debug("No valid move: Switch would violate fixed ring size");
        return false;
      }
    }
  }

  // Check move will not violate ring size limits
  if (networkB.nodes[ringNode1].numConnections() == minRingSize ||
      networkB.nodes[ringNode2].numConnections() == minRingSize ||
      networkB.nodes[ringNode3].numConnections() == maxRingSize ||
      networkB.nodes[ringNode4].numConnections() == maxRingSize) {
    logger->debug("No valid move: Switch would violate ring size limits");
    return false;
  }
  return true;
}

/**
 * @brief Find a common base node connection between a base node and ring node
 * excluding a given node
 * @param baseNode ID of the base node
 * @param ringNode ID of the ring node
 * @param excludeNode ID of the base node to be excluded
 * @return ID of the common connection
 * @throw std::runtime_error if the associated node cannot be found
 */
uint16_t LinkedNetwork::findCommonConnection(const uint16_t baseNode,
                                             const uint16_t ringNode,
                                             const uint16_t excludeNode) const {
  // Find node that shares baseNode and ringNode but is not excludeNode
  std::unordered_set<uint16_t> commonConnections =
      intersectSets(networkA.nodes[baseNode].netConnections,
                    networkB.nodes[ringNode].dualConnections);
  commonConnections.erase(excludeNode);
  if (commonConnections.empty()) {
    throw std::runtime_error(
        std::format("Could not find common base node for base node {} and ring "
                    "node {} excluding node {}",
                    baseNode, ringNode, excludeNode));
  }
  return *commonConnections.begin();
}
/**
 * @brief Find a common ring connection between two base nodes that exlcudes a
 * given node
 * @param baseNode1 ID of the first base node
 * @param baseNode2 ID of the second base node
 * @param excludeNode ID of the ring node to be excluded
 * @return ID of the common ring connection
 * @throw std::runtime_error if the associated node cannot be found
 */
uint16_t LinkedNetwork::findCommonRing(const uint16_t baseNode1,
                                       const uint16_t baseNode2,
                                       const uint16_t excludeNode) const {
  // Find node that shares baseNode1 and baseNode2 but is not excludeNode
  std::unordered_set<uint16_t> commonRings =
      intersectSets(networkA.nodes[baseNode1].dualConnections,
                    networkA.nodes[baseNode2].dualConnections);
  commonRings.erase(excludeNode);
  if (commonRings.empty()) {
    std::string dualConnectionsString1 =
        containerToString(networkA.nodes[baseNode1].dualConnections);
    std::string dualConnectionsString2 =
        containerToString(networkA.nodes[baseNode2].dualConnections);
    throw std::runtime_error(std::format(
        "Could not find common ring node for base node {} and base "
        "node {} excluding ring node {}. Sets: {}\t{} Intersection: {}",
        baseNode1, baseNode2, excludeNode, dualConnectionsString1,
        dualConnectionsString2, containerToString(commonRings)));
  }
  return *commonRings.begin();
}

/**
 * @brief Switch the BSS network by breaking and making bonds
 * @param bondBreaks the bonds to break (vector of pairs)
 * @param bondMakes the bonds to make (vector of pairs)
 */
void LinkedNetwork::applyMove(
    const std::array<std::array<uint16_t, 2>, 2> &bondBreaks,
    const std::array<std::array<uint16_t, 2>, 2> &ringBondBreakMake) {
  uint16_t atom1 = bondBreaks[0][0];
  uint16_t atom2 = bondBreaks[1][0];
  uint16_t atom4 = bondBreaks[1][1];
  uint16_t atom5 = bondBreaks[0][1];

  uint16_t ringNode1 = ringBondBreakMake[0][0];
  uint16_t ringNode2 = ringBondBreakMake[0][1];
  uint16_t ringNode3 = ringBondBreakMake[1][0];
  uint16_t ringNode4 = ringBondBreakMake[1][1];

  auto &nodeA1 = networkA.nodes[atom1];
  auto &nodeA2 = networkA.nodes[atom2];
  auto &nodeA4 = networkA.nodes[atom4];
  auto &nodeA5 = networkA.nodes[atom5];

  auto &nodeB1 = networkB.nodes[ringNode1];
  auto &nodeB2 = networkB.nodes[ringNode2];
  auto &nodeB3 = networkB.nodes[ringNode3];
  auto &nodeB4 = networkB.nodes[ringNode4];

  // A-A connectivities
  setReplace(nodeA1.netConnections, atom5, atom4);
  setReplace(nodeA2.netConnections, atom4, atom5);
  setReplace(nodeA4.netConnections, atom2, atom1);
  setReplace(nodeA5.netConnections, atom1, atom2);

  // A-B connectvities
  setReplace(nodeA1.dualConnections, ringNode1, ringNode4);
  setReplace(nodeA2.dualConnections, ringNode2, ringNode3);

  // B-B connectivities
  nodeB1.netConnections.erase(ringNode2);
  nodeB2.netConnections.erase(ringNode1);
  nodeB3.netConnections.insert(ringNode4);
  nodeB4.netConnections.insert(ringNode3);

  // B-A connectivities
  nodeB1.dualConnections.erase(atom1);
  nodeB2.dualConnections.erase(atom2);

  nodeB3.dualConnections.insert(atom2);
  nodeB4.dualConnections.insert(atom1);
}

/**
 * @brief Restores the initial state of the network by assigning Node objects to
 * their initial states
 * @param initialInvolvedNodesA the initial state of the nodes in network A
 * @param initialInvolvedNodesB the initial state of the nodes in network B
 */
void LinkedNetwork::revertMove(const std::vector<Node> &initialInvolvedNodesA,
                               const std::vector<Node> &initialInvolvedNodesB) {
  // Revert changes to descriptors due to breaking connections
  for (const Node &node : initialInvolvedNodesA) {
    networkA.nodes[node.id] = node;
  }
  for (const Node &node : initialInvolvedNodesB) {
    networkB.nodes[node.id] = node;
  }
}

/**
 * @brief Sets the coordinates of nodes in the network. Ring network coordinates
 * are recalculated.
 * @param coords Coordinates of the base nodes
 * @throw std::invalid_argument if the number of coordinates does not match the
 * number of nodes in network A
 */
void LinkedNetwork::pushCoords(
    const std::vector<std::array<double, 2>> &coords) {
  if (coords.size() != networkA.nodes.size()) {
    throw std::invalid_argument(
        "Number of coordinates does not match number of nodes in network A");
  }

  // Sync A coordinates
  for (int i = 0; i < networkA.nodes.size(); ++i) {
    networkA.nodes[i].coord = coords[i];
  }
  // Centre all the rings relative to their dual connections
  networkB.centerByDual(networkA);
}

/**
 * @brief Checks if network is consistent
 * @return true if all connectivities are reciprocated and all base nodes have
 * clockwise neighbours and ring sizes that are not neighbours of fixedRings are
 * within range, false otherwise
 */
bool LinkedNetwork::checkConsistency() {
  logger->info("Checking consistency...");
  logger->info("Checking base-base connections reciprocated...");
  bool consistent = this->networkA.checkConnectionsReciprocated(this->logger);
  logger->info("Checking ring-ring connections reciprocated...");
  consistent = this->networkB.checkConnectionsReciprocated(this->logger);
  logger->info("Checking base-ring connections reciprocated...");
  consistent =
      this->networkA.checkConnectionsReciprocated(this->networkB, this->logger);
  logger->info("Checking ring-base connections reciprocated...");
  consistent =
      this->networkB.checkConnectionsReciprocated(this->networkA, this->logger);
  logger->info("Checking ring sizes...");
  consistent = this->networkB.checkDegreeLimits(
      this->minRingSize, this->maxRingSize, this->logger);
  return consistent;
}

/**
 * @brief Wraps coordinates out of bounds back into the periodic box, only for 2
 * dimensions
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
 * @brief Checks if the bonds to a given node are within the maximum bond length
 * @param nodeID ID of the node to check
 * @param coords Coordinates of all nodes as a 1D vector of coordinate pairs
 * @return true if all bonds are within the maximum bond length, false otherwise
 */
bool LinkedNetwork::checkBondLengths(
    const uint16_t nodeID,
    const std::vector<std::array<double, 2>> &coords) const {
  for (const auto &neighbourID : networkA.nodes[nodeID].netConnections) {
    std::array<double, 2> pbcVec =
        pbcArray(coords[nodeID], coords[neighbourID], dimensions);
    if (std::sqrt(pbcVec[0] * pbcVec[0] + pbcVec[1] * pbcVec[1]) >
        maximumBondLength) {
      return false;
    }
  }
  return true;
}

/**
 * @brief Checks if the bonds to the given nodes are within the maximum bond
 * length
 * @param nodeIDs IDs of the nodes to check
 * @param coords Coordinates of all nodes as a 1D vector of coordinate pairs
 * @return true if all bonds are within the maximum bond length, false otherwise
 */
bool LinkedNetwork::checkBondLengths(
    const std::unordered_set<uint16_t> &nodeIDs,
    const std::vector<std::array<double, 2>> &coords) const {
  return std::ranges::all_of(nodeIDs, [this, &coords](int nodeID) {
    return checkBondLengths(nodeID, coords);
  });
}

/**
 * @brief Gets the order of a given array of 4 ring IDs using the centre of all
 * of their coordinates
 * @param ringNodeIDs IDs of the rings to get the order of
 * @return clockwise or anti-clockwise
 */
Direction LinkedNetwork::getRingsDirection(
    const std::array<uint16_t, 4> &ringNodeIDs) const {
  const std::array<double, 2> midCoord = this->networkB.getAverageCoordsPBC(
      {ringNodeIDs[0], ringNodeIDs[1], ringNodeIDs[2], ringNodeIDs[3]});
  uint8_t timesDecreased = 0;
  double prevAngle = getClockwiseAnglePBC(
      midCoord, networkB.nodes[ringNodeIDs.back()].coord, dimensions);
  for (uint16_t ringNodeID : ringNodeIDs) {
    double angle = getClockwiseAnglePBC(
        midCoord, this->networkB.nodes[ringNodeID].coord, dimensions);
    if (angle < prevAngle) {
      timesDecreased++;
      if (timesDecreased == 2) {
        return Direction::ANTICLOCKWISE;
      }
    }
    prevAngle = angle;
  }
  return Direction::CLOCKWISE;
}

/**
 * @brief Gets the areas of all rings in the network
 * @return Vector of areas of all rings
 */
std::vector<double> LinkedNetwork::getRingAreas() const {
  std::vector<std::array<double, 2>> currentCoords = this->networkA.getCoords();
  std::vector<double> ringAreas(networkB.nodes.size());

  std::ranges::transform(
      networkB.nodes, ringAreas.begin(),
      [this, &currentCoords](const Node &ringNode) {
        // For each ring node, get the relative coordinates of the connected
        // base nodes
        std::vector<std::array<double, 2>> baseNodeCoords(
            ringNode.dualConnections.size());
        std::ranges::transform(
            ringNode.dualConnections, std::back_inserter(baseNodeCoords),
            [this, &currentCoords, &ringNode](int baseNodeID) {
              return pbcArray(ringNode.coord, currentCoords[baseNodeID],
                              dimensions);
            });
        // Calculate the polygon area of these relative coordinates
        return calculatePolygonArea(baseNodeCoords);
      });
  return ringAreas;
}

bool LinkedNetwork::checkAngles(
    const std::unordered_set<uint16_t> &nodeIDs,
    const std::vector<std::array<double, 2>> &potentialCoords) const {
  for (const uint16_t nodeID : nodeIDs) {
    const Node &node = this->networkA.nodes[nodeID];
    // Construct connection coordinates
    std::vector<std::array<double, 2>> connectionCoords(
        node.netConnections.size());
    for (const uint16_t connectionID : node.netConnections) {
      connectionCoords.push_back(potentialCoords[connectionID]);
    }
    // Check angles with periodic boundary conditions
    if (checkAnglesPBC(node.coord, connectionCoords, this->dimensions, 0,
                       this->maximumAngle)) {
      continue;
    }
    // Check failed
    logger->debug("Angle check failed for node {} with connections: {}", nodeID,
                  containerToString(node.netConnections));
    return false;
  }
  return true;
}