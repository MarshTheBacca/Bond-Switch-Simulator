#ifndef NL_NODE_H
#define NL_NODE_H

#include <array>
#include <cmath>
#include <set>
#include <string>

struct Node {
  // Unique ID number
  int id = 0;
  // Location of the node in 2D euclidean space
  std::array<double, 2> coord{0.0, 0.0};
  // IDs of connected nodes in the same network
  std::set<int> netConnections{};
  // IDs of connected nodes in the dual network
  std::set<int> dualConnections{};

  /**
   * @brief Construct a node with an ID of 0, no connections and coordinate at
   * (0, 0)
   */
  Node() = default;

  /**
   * @brief Construct node with a given ID with no connections and coordinate at
   * (0, 0)
   * @param nodeId The ID of the node
   */
  explicit Node(const int nodeID);

  /**
   * @brief Constructor with a given node ID and coordinates and no connections
   * @param nodeID The ID of the node
   * @param coord The coordinates of the node
   */
  Node(const int nodeID, const std::array<double, 2> &coord);

  /**
   * @brief Constructor with a given node ID, coordinates and connections
   * @param nodeID The ID of the node
   * @param coord The coordinates of the node
   * @param netConnections The connections to nodes in the network
   * @param dualConnections The connections to nodes in the dual network
   */
  Node(const int nodeID, const std::array<double, 2> &coord,
       const std::set<int> &netConnections,
       const std::set<int> &dualConnections);

  /**
   * @brief Convert the node to a string
   */
  std::string toString() const;

  /**
   * @brief Calculate the distance between this node and a given coordinate
   * @param coordinate The coordinate to calculate the distance to
   * @return distance
   */
  inline double distanceFrom(const std::array<double, 2> &coordinate) {
    return sqrt(pow(coord[0] - coordinate[0], 2) +
                pow(coord[1] - coordinate[1], 2));
  }

  /**
   * @brief Get the node's number of connections
   * @return Number of connections
   */
  inline int numConnections() const {
    return static_cast<int>(netConnections.size());
  }

  /**
   * @brief Get the node's number of dual connections
   * @return Number of dual connections
   */
  inline int numDualConnections() const {
    return static_cast<int>(netConnections.size());
  }
};

#endif // NL_NODE_H
