#include "node.h"

// Default constructor
Node::Node() {
    id = -1;
    crd = VecF<double>(1);
    netCnxs = VecR<int>(1);
    dualCnxs = VecR<int>(1);
    auxCnxs = VecR<int>(1);
}

// Construct with maximum number of connections, initialise with 0
Node::Node(int nodeId, int maxNetCnxs, int maxDualCnxs, int maxAuxCnxs) : id(nodeId), netCnxs(0, maxNetCnxs),
                                                                          dualCnxs(0, maxDualCnxs), auxCnxs(0, maxAuxCnxs) {
    oxyCnxs = VecR<int>(8);
}

std::ostream &operator<<(std::ostream &os, Node &node) {
    // Replace this with your actual implementation
    os << "Node " << node.id << "at " << node.crd[0] << ", " << node.crd[1] << " with neighbours: ";
    for (int i = 0; i < node.netCnxs.n; ++i)
        os << node.netCnxs[i] << " ";
    os << " and ring neighbours: ";
    for (int i = 0; i < node.dualCnxs.n; ++i)
        os << node.dualCnxs[i] << " ";
    return os;
}

std::ostream &operator<<(std::ostream &os, const Node &node) {
    // Replace this with your actual implementation
    os << "Node " << node.id << "at " << node.crd[0] << ", " << node.crd[1] << " with neighbours: ";
    for (int i = 0; i < node.netCnxs.n; ++i)
        os << node.netCnxs[i] << " ";
    os << " and ring neighbours: ";
    for (int i = 0; i < node.dualCnxs.n; ++i)
        os << node.dualCnxs[i] << " ";
    return os;
}

std::string Node::toString() {
    std::string str = "Node " + std::to_string(id) + " at " + std::to_string(crd[0]) + ", " + std::to_string(crd[1]) +
                      " with neighbours: ";
    for (int i = 0; i < netCnxs.n; ++i)
        str += std::to_string(netCnxs[i]) + " ";
    str += " and ring neighbours: ";
    for (int i = 0; i < dualCnxs.n; ++i)
        str += std::to_string(dualCnxs[i]) + " ";
    return str;
}