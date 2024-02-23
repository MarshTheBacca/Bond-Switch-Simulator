#include "node.h"

// Default constructor
Node::Node() : id(-1), crd(1), netCnxs(1), dualCnxs(1) {
}

// Construct with maximum number of connections, initialise with 0
Node::Node(const int &nodeId, const int &maxNetCnxs, const int &maxDualCnxs) : id(nodeId), netCnxs(0, maxNetCnxs),
                                                                               dualCnxs(0, maxDualCnxs) {
}

std::ostream &operator<<(std::ostream &os, const Node &node) {
    // Replace this with your actual implementation
    os << "Node " << node.id << "at " << node.crd[0] << ", " << node.crd[1] << " with neighbours: ";
    std::for_each(node.netCnxs.begin(), node.netCnxs.end(), [&os](int i) { os << i << " "; });
    os << " and ring neighbours: ";
    std::for_each(node.dualCnxs.begin(), node.dualCnxs.end(), [&os](int i) { os << i << " "; });
    return os;
}

std::string Node::toString() {
    std::string str = "Node " + std::to_string(id) + " at " + std::to_string(crd[0]) + ", " + std::to_string(crd[1]) +
                      " with neighbours: ";
    std::for_each(netCnxs.begin(), netCnxs.end(), [&str](int i) { str += std::to_string(i) + " "; });
    str += " and ring neighbours: ";
    std::for_each(dualCnxs.begin(), dualCnxs.end(), [&str](int i) { str += std::to_string(i) + " "; });
    return str;
}
