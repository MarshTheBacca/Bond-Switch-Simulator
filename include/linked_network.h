// Linked reciprocal networks - network and dual pair

#ifndef NL_LINKED_NETWORK_H
#define NL_LINKED_NETWORK_H
#include "input_data.h"
#include "lammps_object.h"
#include "metropolis.h"
#include "network.h"
#include <array>
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
  std::unordered_map<int, int> fixedRings;
  std::unordered_set<int> fixedNodes; // IDs of the fixed nodes

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
  std::tuple<int, int, int, int> pickRandomConnection();
  int assignValues(int randNodeCoordination,
                   int randNodeConnectionCoordination) const;

  int findCommonConnection(const int baseNode, const int ringNode,
                           const int excludeNode) const;
  int findCommonRing(const int baseNode1, const int baseNode2,
                     const int excludeNode) const;

  void monteCarloSwitchMoveLAMMPS(const double &temperature);
  void rejectMove(const std::vector<Node> &initialInvolvedNodesA,
                  const std::vector<Node> &initialInvolvedNodesB,
                  const std::vector<int> &bondBreaks,
                  const std::vector<int> &bondMakes,
                  const std::vector<int> &angleBreaks,
                  const std::vector<int> &angleMakes);

  bool checkConsistency();

  void write() const;

  void pushCoords(const std::vector<double> &coords);
  void showCoords(const std::vector<double> &coords) const;
  void wrapCoords(std::vector<double> &coords) const;

  bool genSwitchOperations(int baseNode1, int baseNode2, int ringNode1,
                           int ringNode2, std::vector<int> &bondBreaks,
                           std::vector<int> &bondMakes,
                           std::vector<int> &angleBreaks,
                           std::vector<int> &angleMakes,
                           std::vector<int> &ringBondBreakMake,
                           std::unordered_set<int> &convexCheckIDs);

  void switchNetMCGraphene(const std::vector<int> &bondBreaks,
                           const std::vector<int> &ringBondBreakMake);
  void revertNetMCGraphene(const std::vector<Node> &initialInvolvedNodesA,
                           const std::vector<Node> &initialInvolvedNodesB);

  std::tuple<std::vector<double>, std::vector<double>>
  rotateBond(const int atomID1, const int atomID2,
             const Direction &direct) const;
  Direction getRingsDirection(const std::vector<int> &ringNodeIDs) const;

  bool checkClockwiseNeighbours(const int nodeID) const;
  bool checkClockwiseNeighbours(const int nodeID,
                                const std::vector<double> &coords) const;
  bool checkAllClockwiseNeighbours() const;
  void arrangeNeighboursClockwise(const int nodeID,
                                  const std::vector<double> &coords);
  void arrangeNeighboursClockwise(const std::unordered_set<int> &nodeIDs,
                                  const std::vector<double> &coords);

  bool checkAnglesWithinRange(const std::vector<double> &coords);
  bool checkAnglesWithinRange(const std::unordered_set<int> &nodeIDs,
                              const std::vector<double> &coords);
  bool checkBondLengths(const int nodeID,
                        const std::vector<double> &coords) const;
  bool checkBondLengths(const std::unordered_set<int> &nodeIDs,
                        const std::vector<double> &coords) const;

  std::map<int, double> getRingSizes() const;
  std::vector<double> getRingAreas() const;
};

#endif // NL_LINKED_NETWORK_H
