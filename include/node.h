// Node in network containing coordinate, connections to nodes in network and connections to nodes in dual

#ifndef NL_NODE_H
#define NL_NODE_H

#include "vec_func.h"
#include "vecf.h"
#include "vecr.h"
#include <iostream>

// Node in network, aka vertex in graph
class Node {
  public:
    // Data members
    int id;             // unique id used by other nodes for connections
    VecF<double> crd;   // coordinate
    VecR<int> netCnxs;  // connections to nodes in network
    VecR<int> dualCnxs; // connections to nodes in dual
    VecR<int> auxCnxs;  // free form additional connections to be repurposed as needed
    VecR<int> oxyCnxs;  // hack as cannot return arrays in c++11

    // Constructors
    Node();
    Node(int nodeId, int maxNetCnxs, int maxDualCnxs, int maxAuxCnxs);
    inline double distanceFrom(VecF<double> &coordinate) {
        return sqrt(pow(crd[0] - coordinate[0], 2) + pow(crd[1] - coordinate[1], 2));
    }

    // Destructor
    ~Node() = default;

    // Member functions

    std::string toString();
};

std::ostream &operator<<(std::ostream &os, Node &node);
std::ostream &operator<<(std::ostream &os, const Node &node);

#endif // NL_NODE_H
