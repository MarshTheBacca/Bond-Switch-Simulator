#include "node.h"
#include "random_number_generator.h"
#include <algorithm>
#include <array>
#include <ranges>
#include <sstream>
#include <string>
#include <unordered_set>

/**
 * @brief Construct node with a given ID with no connections and coordinate at
 * (0, 0)
 * @param nodeId The ID of the node
 */
Node::Node(const uint16_t nodeId) : id(nodeId) {}

/**
 * @brief Constructor with a given node ID and coordinates and no connections
 * @param nodeID The ID of the node
 * @param coord The coordinates of the node
 */
Node::Node(const uint16_t nodeID, const std::array<double, 2> &coord)
    : id(nodeID), coord(coord) {}

/**
 * @brief Constructor with a given node ID, coordinates and connections
 * @param nodeID The ID of the node
 * @param coord The coordinates of the node
 * @param netConnections The connections to nodes in the network
 * @param dualConnections The connections to nodes in the dual network
 */
Node::Node(const uint16_t nodeID, const std::array<double, 2> &coord,
           const std::unordered_set<uint16_t> &netConnections,
           const std::unordered_set<uint16_t> &dualConnections)
    : id(nodeID), coord(coord), netConnections(netConnections),
      dualConnections(dualConnections) {}

/**
 * @brief Convert the node to a string
 */
std::string Node::toString() const {
  std::ostringstream oss;
  oss << std::format("Node {} at {}, {} with neighbours: ", id, coord[0],
                     coord[1]);
  std::ranges::for_each(netConnections,
                        [&oss](int i) { oss << std::format("{} ", i); });
  oss << " and ring neighbours: ";
  std::ranges::for_each(dualConnections,
                        [&oss](int i) { oss << std::format("{} ", i); });
  return oss.str();
}

uint16_t Node::getRandomNetConnectionID() const {
  return *RandomNumberGenerator::getInstance().getRandomElement(
      this->netConnections.begin(), this->netConnections.end());
}

uint16_t Node::getRandomDualConnectionID() const {
  return *RandomNumberGenerator::getInstance().getRandomElement(
      this->dualConnections.begin(), this->dualConnections.end());
}

Node &Node::operator=(const Node &other) {
  if (this != &other) {
    this->coord = other.coord;
    this->netConnections = other.netConnections;
    this->dualConnections = other.dualConnections;
  }
  return *this;
}