#ifndef LAMMPS_OBJECT_H
#define LAMMPS_OBJECT_H

#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>

#include "vec_func.h"
#include "vecf.h"
#include "vecr.h"

#include "mpi.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"

#include "network.h"

#include "atom.h"
#include "input.h"
#include "lammps.h" // these are LAMMPS include files
#include "library.h"

#include <spdlog/spdlog.h>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

class LammpsObject {
  public:
    Network networkA;
    Network networkB;

    void *handle;
    int version;
    int natoms = 0;
    int nbonds = 0;
    int nangles = 0;
    double *bonds = nullptr;

    std::string prefixFolderIn;
    std::string prefixFolderOut;
    std::string prefixOut;
    std::vector<int> angleHelper;

    LammpsObject();
    LammpsObject(const std::string &selector, const std::string &inputFolder, const LoggerPtr logger);

    void write_data(const std::string &structureName);
    void write_restart(const std::string &structureName);
    void finaliseLAMMPSObject(const std::string &structureName);

    std::vector<double> getCoords(int dim);
    void getCoords(std::vector<double> &coords, int dim);
    void setCoords(std::vector<double> &newCoords, int dim);

    void breakBond(int atom1, int atom2, int type, LoggerPtr logger);
    void formBond(int atom1, int atom2, int type, LoggerPtr logger);
    void breakAngle(int atom1, int atom2, int atom3, LoggerPtr logger);
    void formAngle(int atom1, int atom2, int atom3, LoggerPtr logger);
    void breakAndFormAngles(int atom3, int atom1, int atom5, int atom4, int e1, int e11, int d1, int d11);
    std::pair<int, int> identifyAtoms(int atomA, int atomB, Network networkAArg, LoggerPtr logger);

    void switchGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                        const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes, LoggerPtr logger);
    void revertGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                        const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes, LoggerPtr logger);
    std::vector<int> getAngles();
    void showAngles(const int &numLines, LoggerPtr logger);
    bool checkAngleUnique(const int &atom1, const int &atom2, const int &atom3);

    void minimiseNetwork();
    double getPotentialEnergy();
};

#endif // LAMMPS_OBJECT_H
