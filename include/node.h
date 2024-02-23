// Node in network containing coordinate, connections to nodes in network and connections to nodes in dual

#ifndef NL_NODE_H
#define NL_NODE_H

#include "vector_tools.h"
#include <iostream>

// Node in network, aka vertex in graph
class Node {
  public:
    // Data members
    int id;                    // unique id used by other nodes for connections
    std::vector<double> crd;   // coordinate
    std::vector<int> netCnxs;  // connections to nodes in network
    std::vector<int> dualCnxs; // connections to nodes in dual

    // Constructors
    Node();
    Node(const int &nodeID, const int &maxNetCnxs, const int &maxDualCnxs);

    inline double distanceFrom(const std::vector<double> &coordinate) {
        return sqrt(pow(crd[0] - coordinate[0], 2) + pow(crd[1] - coordinate[1], 2));
    }

    // Member functions
    std::string toString();
};

std::ostream &operator<<(std::ostream &os, const Node &node);

#endif // NL_NODE_H
