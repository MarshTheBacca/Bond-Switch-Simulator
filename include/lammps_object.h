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
public:
    Network networkA;
    Network networkB;
    Network networkT;

    void *handle;
    int version;
    int natoms = 0;
    int nbonds = 0;
    int nangles = 0;
    double *coords = nullptr;
    double *bonds = nullptr;
    double *angles = nullptr;

    std::string prefixFolderIn;
    std::string prefixFolderOut;
    std::string prefixOut;

    LammpsObject();
    LammpsObject(const std::string &selector, const std::string &inputFolder, const LoggerPtr logger);

    void write_data(const std::string &structureName);
    void write_restart(const std::string &structureName);
    void finaliseLAMMPSObject(const std::string &structureName);

    double pbx();
    double pby();
    double pbz();

    double *fetchCrds(int dim);
    void pushCrds(int dim, double *old_coords);

    double *fetchBonds();

    void breakBond(int atom1, int atom2, int type, LoggerPtr logger);
    void formBond(int atom1, int atom2, int type, LoggerPtr logger);
    void breakAngle(int atom1, int atom2, int atom3);
    void formAngle(int atom1, int atom2, int atom3);

    void switchGraphene(VecF<int> switchIdsA, Network networkA, LoggerPtr logger);
    void revertGraphene(VecF<int> switchIdsA, Network networkA, LoggerPtr logger);

    void switchTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger);
    void revertTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger);

    void switchBilayer(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger);
    void revertBilayer(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger);

    void switchBonds(VecF<int> switchIdsA, VecF<int> switchIdsT);
    void revertBonds(VecF<int> switchIdsA, VecF<int> switchIdsT);

    VecF<int> GlobalPotentialMinimisation();
    double globalPotentialEnergy();

    double probeE();
};

#endif // LAMMPS_OBJECT_H
