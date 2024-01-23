//
// Created by olwhi on 24/07/2023.
//

#ifndef NETMC_LAMMPS_C_INTERFACE_H
#define NETMC_LAMMPS_C_INTERFACE_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <random>

#include "vecf.h"
#include "vecr.h"
#include "vec_func.h"

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "mpi.h"

#include "network.h"
// #include "linked_network.h"
// #include "node.h"

#include "lammps.h" // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

using namespace std;

class LammpsObject
{

private:
public:
    LammpsObject();
    LammpsObject(int selector, string prefixin, string prefixout);

    Network networkA, networkB, networkT;

    void *handle;
    int version;
    int natoms = 0, nbonds = 0, nangles = 0;
    double *coords = NULL;
    double *bonds = NULL;
    double *angles = NULL;

    string prefixFolderIn, prefixFolderOut;

    bool globalVerbose = false;
    string prefixOut;

    int initialiseLammpsObject();
    int write_data(int selector);
    int write_restart(int selector);
    int finaliseLammpsObject(int selector);

    void runInput(string fname);
    void getatominfo(int dim);

    double pbx();
    double pby();
    double pbz();

    double *fetchCrds(int dim);
    void pushCrds(int dim, double *old_coords);

    double *fetchBonds();

    int getnAtoms();

    void breakBond(int atom1, int atom2, int type);
    void formBond(int atom1, int atom2, int type);
    void breakAngle(int atom1, int atom2, int atom3);
    void formAngle(int atom1, int atom2, int atom3);

    void switchGraphene(VecF<int> switchIdsA, Network networkA);
    void revertGraphene(VecF<int> switchIdsA, Network networkA);

    void switchTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT, Network networkT);
    void revertTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);

    void switchBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);
    void revertBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);

    void switchBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);
    void revertBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);

    VecF<int> GlobalPotentialMinimisation();
    double GlobalPotentialEnergy();

    double probeE();
};

#endif // NETMC_LAMMPS_C_INTERFACE_H
