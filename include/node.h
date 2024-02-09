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

    // Constructors, copy constructor, assignment operator
    Node();
    Node(int nodeId, int maxNetCnxs, int maxDualCnxs, int maxAuxCnxs);
    Node(const Node &source);
    Node &operator=(const Node &source);
    double distanceFrom(VecF<double> &crd2);
};

#endif // NL_NODE_H
