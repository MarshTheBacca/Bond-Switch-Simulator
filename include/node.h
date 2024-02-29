#ifndef NL_NODE_H
#define NL_NODE_H

#include "vector_tools.h"
#include <iostream>

struct Node {
    // Member variables
    int id = 0;                        // unique id used by other nodes for connections
    std::vector<double> crd{0.0, 0.0}; // coordinate
    std::vector<int> netCnxs;          // connections to nodes in network
    std::vector<int> dualCnxs;         // connections to nodes in dual network

    // Constructors
    Node();
    explicit Node(const int &nodeID);
    Node(const int &nodeID, const std::vector<double> &crd);
    Node(const int &nodeID, const std::vector<double> &crd, const std::vector<int> &netCnxs, const std::vector<int> &dualCnxs);

    // Methods
    std::string toString();

    /**
     * @brief Calculate the distance between this node and a given coordinate
     * @param coordinate The coordinate to calculate the distance to
     * @return distance
     */
    inline double distanceFrom(const std::vector<double> &coordinate) {
        return sqrt(pow(crd[0] - coordinate[0], 2) + pow(crd[1] - coordinate[1], 2));
    }
};

#endif // NL_NODE_H
