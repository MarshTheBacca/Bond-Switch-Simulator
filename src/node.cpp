#include "node.h"
#include <algorithm>
#include <array>
#include <ranges>
#include <set>
#include <sstream>
#include <string>

Node::Node(const int nodeId) : id(nodeId) {}

Node::Node(const int nodeID, const std::array<double, 2> &coord)
    : id(nodeID), coord(coord) {}

Node::Node(const int nodeID, const std::array<double, 2> &coord,
           const std::set<int> &netConnections,
           const std::set<int> &dualConnections)
    : id(nodeID), coord(coord), netConnections(netConnections),
      dualConnections(dualConnections) {}

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
