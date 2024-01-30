#ifndef LAMMPS_OBJECT_H
#define LAMMPS_OBJECT_H

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

#include "lammps.h" // these are LAMMPS include files
#include "input.h"
#include "atom.h"
#include "library.h"

#include <spdlog/spdlog.h>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

class LammpsObject
{

private:
public:
    Network networkA, networkB, networkT;

    void *handle;
    int version;
    int natoms = 0, nbonds = 0, nangles = 0;
    double *coords = NULL;
    double *bonds = NULL;
    double *angles = NULL;

    std::string prefixFolderIn, prefixFolderOut;
    std::string prefixOut;

    LammpsObject();
    LammpsObject(std::string selector, std::string inputFolder, LoggerPtr logger);

    int write_data(int selector);
    int write_restart(int selector);
    int finaliseLammpsObject(int selector);

    void runInput(std::string fname);
    void getatominfo(int dim);

    double pbx();
    double pby();
    double pbz();

    double *fetchCrds(int dim);
    void pushCrds(int dim, double *old_coords);

    double *fetchBonds();

    int getnAtoms();

    void breakBond(int atom1, int atom2, int type, LoggerPtr logger);
    void formBond(int atom1, int atom2, int type, LoggerPtr logger);
    void breakAngle(int atom1, int atom2, int atom3);
    void formAngle(int atom1, int atom2, int atom3);

    void switchGraphene(VecF<int> switchIdsA, Network networkA, LoggerPtr logger);
    void revertGraphene(VecF<int> switchIdsA, Network networkA, LoggerPtr logger);

    void switchTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT, Network networkT, LoggerPtr logger);
    void revertTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT, LoggerPtr logger);

    void switchBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT, LoggerPtr logger);
    void revertBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT, LoggerPtr logger);

    void switchBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);
    void revertBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT);

    VecF<int> GlobalPotentialMinimisation();
    double GlobalPotentialEnergy();

    double probeE();
};

#endif // LAMMPS_OBJECT_H
