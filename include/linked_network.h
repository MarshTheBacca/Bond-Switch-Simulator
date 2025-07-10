#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H

#include "input_data.h"
#include "lammps_manager.h"
#include "network.h"
#include "stats.h"
#include "switch_move.h"
#include "types.h"
#include <array>
#include <cstdint>
#include <optional>
#include <spdlog/spdlog.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

struct LinkedNetwork {
  // Data members
  Network networkB; // Ring network
  int minRingSize;  // Minimum coordination number of ring network
  int maxRingSize;  // Maximum coordination number of ring network

  Network networkA; // Base network

  std::array<double, 2> dimensions; // Periodic boundary of network, xlo = ylo =
                                    // 0, so dimensions = [xhi, yhi]
  std::array<double, 2> centreCoords; // Centre of network = [xhi / 2, yhi / 2]

  LAMMPSManager lammpsManager; // LAMMPSManager for network
  double energy;               // The current energy of the system

  SelectionType selectionType; // Either 'weighted' or 'random'
  double weightedDecay;        // decay factor for weighted monte carlo
  double maximumBondLength;    // Maximum bond length
  double maximumAngle;         // Maximum angle between atoms
  bool writeMovie;             // Write movie file or not

  // Map of Fixed Ring IDs to their size
  std::unordered_map<uint16_t, size_t> fixedRings;
  std::unordered_set<uint16_t> fixedNodes; // IDs of the fixed nodes

  // Statistics
  Stats stats;

  LoggerPtr logger; // Logger
  std::vector<double> weights;

  // Constructors
  /**
   * @brief Default constructor
   */
  LinkedNetwork() = default;
  LinkedNetwork(const InputData &inputData, const LoggerPtr &logger);

  void findFixedRings(const std::string &flePath);
  void findFixedNodes();

  void rescale(double scaleFactor);
  void updateWeights();

  std::array<std::array<uint16_t, 2>, 2> pickRandomConnection();

  uint16_t findCommonConnection(const uint16_t baseNode,
                                const uint16_t ringNode,
                                const uint16_t excludeNode) const;
  uint16_t findCommonRing(const uint16_t baseNode1, const uint16_t baseNode2,
                          const uint16_t excludeNode) const;

  void performBondSwitch(const double temperature);
  void acceptMove(const std::vector<std::array<double, 2>> &newCoords,
                  const double initialEnergy, const double finalEnergy);
  void rejectMove(const std::vector<Node> &initialInvolvedNodesA,
                  const std::vector<Node> &initialInvolvedNodesB,
                  const SwitchMove &switchMove);

  bool checkConsistency();

  void write() const;

  void pushCoords(const std::vector<std::array<double, 2>> &coords);
  void showCoords(const std::vector<std::array<double, 2>> &coords) const;
  void wrapCoords(std::vector<double> &coords) const;

  void
  applyMove(const std::array<std::array<uint16_t, 2>, 2> &bondBreaks,
            const std::array<std::array<uint16_t, 2>, 2> &ringBondBreakMake);
  void revertMove(const std::vector<Node> &initialInvolvedNodesA,
                  const std::vector<Node> &initialInvolvedNodesB);

  Direction getRingsDirection(const std::array<uint16_t, 4> &ringNodeIDs) const;

  bool checkBondLengths(const uint16_t nodeID,
                        const std::vector<std::array<double, 2>> &coords) const;
  bool checkBondLengths(const std::unordered_set<uint16_t> &nodeIDs,
                        const std::vector<std::array<double, 2>> &coords) const;

  std::vector<double> getRingAreas() const;

  bool checkConnectivityLimits(const uint16_t ringNode1,
                               const uint16_t ringNode2,
                               const uint16_t ringNode3,
                               const uint16_t ringNode4) const;

  std::optional<SwitchMove> getSwitchMove();
  SwitchMove genSwitchMove(const uint16_t baseNode1, const uint16_t baseNode2,
                           const uint16_t ringNode1,
                           const uint16_t ringNode2) const;
  bool
  checkAngles(const std::unordered_set<uint16_t> &nodeIDs,
              const std::vector<std::array<double, 2>> &potentialCoords) const;
};

#endif // NL_LINKED_NETWORK_H
