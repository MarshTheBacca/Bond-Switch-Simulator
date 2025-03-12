#ifndef NL_NODE_H
#define NL_NODE_H

#include <array>
#include <cmath>
#include <cstdint>
#include <string>
#include <unordered_set>

struct Node {
  // Unique ID number
  const uint16_t id = 0;
  // Location of the node in 2D euclidean space
  std::array<double, 2> coord{0.0, 0.0};
  // IDs of connected nodes in the same network
  std::unordered_set<uint16_t> netConnections{};
  // IDs of connected nodes in the dual network
  std::unordered_set<uint16_t> dualConnections{};

  /**
   * @brief Construct a node with an ID of 0, no connections and coordinate at
   * (0, 0)
   */
  Node() = default;

  explicit Node(const uint16_t nodeID);

  Node(const uint16_t nodeID, const std::array<double, 2> &coord);

  Node(const uint16_t nodeID, const std::array<double, 2> &coord,
       const std::unordered_set<uint16_t> &netConnections,
       const std::unordered_set<uint16_t> &dualConnections);

  Node &operator=(const Node &other);

  std::string toString() const;

  /**
   * @brief Calculate the distance between this node and a given coordinate
   * @param coordinate The coordinate to calculate the distance to
   * @return distance
   */
  inline double distanceFrom(const std::array<double, 2> &coordinate) const {
    return sqrt(pow(coord[0] - coordinate[0], 2) +
                pow(coord[1] - coordinate[1], 2));
  }

  /**
   * @brief Get the node's number of connections
   * @return Number of connections
   */
  inline size_t numConnections() const { return this->netConnections.size(); }

  /**
   * @brief Get the node's number of dual connections
   * @return Number of dual connections
   */
  inline size_t numDualConnections() const {
    return this->dualConnections.size();
  }

  uint16_t getRandomNetConnectionID() const;
  uint16_t getRandomDualConnectionID() const;
};

#endif // NL_NODE_H
