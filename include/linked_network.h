#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H

#include "input_data.h"
#include "lammps_object.h"
#include "metropolis.h"
#include "network.h"
#include "switch_move.h"
#include <array>
#include <cstdint>
#include <optional>
#include <random>
#include <spdlog/spdlog.h>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

enum class Direction { CLOCKWISE, ANTICLOCKWISE };

struct LinkedNetwork {
  // Data members
  Network networkB; // Ring network
  int minRingSize;  // Minimum coordination number of ring network
  int maxRingSize;  // Maximum coordination number of ring network

  Network networkA; // Base network

  std::array<double, 2> dimensions; // Periodic boundary of network, xlo = ylo =
                                    // 0, so dimensions = [xhi, yhi]
  std::array<double, 2> centreCoords; // Centre of network = [xhi / 2, yhi / 2]

  LammpsObject lammpsNetwork; // LAMMPS object for network
  double energy;              // The current energy of the system

  std::vector<double> currentCoords;

  bool isOpenMPIEnabled;          // Whether to use MPI
  SelectionType selectionType;    // Either 'weighted' or 'random'
  std::mt19937 randomNumGen;      // mersenne twister random number generator
  Metropolis metropolisCondition; // monte carlo metropolis condition
  double weightedDecay;           // decay factor for weighted monte carlo
  double maximumBondLength;       // Maximum bond length
  double maximumAngle;            // Maximum angle between atoms
  bool writeMovie;                // Write movie file or not

  // Map of Fixed Ring IDs to their size
  std::unordered_map<uint16_t, size_t> fixedRings;
  std::unordered_set<uint16_t> fixedNodes; // IDs of the fixed nodes

  int numSwitches = 0;            // Number of switches performed
  int numAcceptedSwitches = 0;    // Number of switches accepted
  int failedBondLengthChecks = 0; // Number of failed bond length checks
  int failedAngleChecks = 0;      // Number of failed angle checks
  int failedEnergyChecks = 0;     // Number of failed energy checks

  LoggerPtr logger; // Logger
  std::vector<double> weights;

  // Constructors
  /**
   * @brief Default constructor
   */
  LinkedNetwork() = default;
  LinkedNetwork(const int numRing, const LoggerPtr &logger);
  LinkedNetwork(const InputData &inputData, const LoggerPtr &logger);

  void findFixedRings(const std::string &flePath);
  void findFixedNodes();

  void rescale(double scaleFactor);
  void updateWeights();

  std::array<std::array<uint16_t, 2>, 2> pickRandomConnection();
  int assignValues(int randNodeCoordination,
                   int randNodeConnectionCoordination) const;

  uint16_t findCommonConnection(const uint16_t baseNode,
                                const uint16_t ringNode,
                                const uint16_t excludeNode) const;
  uint16_t findCommonRing(const uint16_t baseNode1, const uint16_t baseNode2,
                          const uint16_t excludeNode) const;

  void performBondSwitch(const double temperature);
  void rejectMove(const std::vector<Node> &initialInvolvedNodesA,
                  const std::vector<Node> &initialInvolvedNodesB,
                  const std::array<std::array<uint16_t, 2>, 2> &bondBreaks,
                  const std::array<std::array<uint16_t, 2>, 2> &bondMakes,
                  const std::array<std::array<uint16_t, 3>, 8> &angleBreaks,
                  const std::array<std::array<uint16_t, 3>, 8> &angleMakes);

  bool checkConsistency();

  void write() const;

  void pushCoords(const std::vector<double> &coords);
  void showCoords(const std::vector<double> &coords) const;
  void wrapCoords(std::vector<double> &coords) const;

  bool
  genSwitchOperations(uint16_t baseNode1, uint16_t baseNode2,
                      uint16_t ringNode1, uint16_t ringNode2,
                      std::array<std::array<uint16_t, 2>, 2> &bondBreaks,
                      std::array<std::array<uint16_t, 2>, 2> &bondMakes,
                      std::array<std::array<uint16_t, 3>, 8> &angleBreaks,
                      std::array<std::array<uint16_t, 3>, 8> &angleMakes,
                      std::array<std::array<uint16_t, 2>, 2> &ringBondBreakMake,
                      std::unordered_set<uint16_t> &convexCheckIDs);

  void
  applyMove(const std::array<std::array<uint16_t, 2>, 2> &bondBreaks,
            const std::array<std::array<uint16_t, 2>, 2> &ringBondBreakMake);
  void revertMove(const std::vector<Node> &initialInvolvedNodesA,
                  const std::vector<Node> &initialInvolvedNodesB);

  std::tuple<std::vector<double>, std::vector<double>>
  rotateBond(const std::array<uint16_t, 2> &bond,
             const Direction &direct) const;
  Direction getRingsDirection(const std::vector<uint16_t> &ringNodeIDs) const;

  bool checkClockwiseNeighbours(const uint16_t nodeID) const;
  bool checkClockwiseNeighbours(const uint16_t nodeID,
                                const std::vector<double> &coords) const;
  bool checkAllClockwiseNeighbours() const;
  void arrangeNeighboursClockwise(const uint16_t nodeID,
                                  const std::vector<double> &coords);
  void arrangeNeighboursClockwise(const std::unordered_set<uint16_t> &nodeIDs,
                                  const std::vector<double> &coords);

  bool checkAnglesWithinRange(const std::vector<double> &coords);
  bool checkAnglesWithinRange(const std::unordered_set<uint16_t> &nodeIDs,
                              const std::vector<double> &coords);
  bool checkBondLengths(const uint16_t nodeID,
                        const std::vector<double> &coords) const;
  bool checkBondLengths(const std::unordered_set<uint16_t> &nodeIDs,
                        const std::vector<double> &coords) const;

  std::map<int, double> getRingSizes() const;
  std::vector<double> getRingAreas() const;

  bool checkConnectivityLimits(const uint16_t ringNode1,
                               const uint16_t ringNode2,
                               const uint16_t ringNode3,
                               const uint16_t ringNode4) const;

  std::optional<SwitchMove> getSwitchMove();
  SwitchMove genSwitchMove(const uint16_t baseNode1, const uint16_t baseNode2,
                           const uint16_t ringNode1,
                           const uint16_t ringNode2) const;
};

#endif // NL_LINKED_NETWORK_H
