#ifndef LAMMPS_OBJECT_H
#define LAMMPS_OBJECT_H

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <sstream>
#include <stdio.h>
#include <unistd.h>


#include "vec_func.h"

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
    LammpsObject(const std::string &selector, const std::string &inputFolder, LoggerPtr logger);

    void write_data(const std::string &structureName);
    void write_restart(const std::string &structureName);

    std::vector<double> getCoords(const int &dim) const;
    void getCoords(std::vector<double> &coords, const int &dim) const;
    void setCoords(std::vector<double> &newCoords, int dim);

    void breakBond(const int &atom1, const int &atom2, const int &type);
    void formBond(const int &atom1, const int &atom2, const int &type);
    void breakAngle(const int &atom1, const int &atom2, const int &atom3, LoggerPtr logger);
    void formAngle(const int &atom1, const int &atom2, const int &atom3, LoggerPtr logger);
    void breakAndFormAngles(int atom3, int atom1, int atom5, int atom4, int e1, int e11, int d1, int d11);
    std::pair<int, int> identifyAtoms(int atomA, int atomB, Network networkAArg, LoggerPtr logger);

    void switchGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                        const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes,
                        const std::vector<double> &rotatedCoord1, const std::vector<double> &rotatedCoord2, LoggerPtr logger);
    void revertGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                        const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes, LoggerPtr logger);
    std::vector<int> getAngles() const;
    void showAngles(const int &numLines, LoggerPtr logger) const;
    bool checkAngleUnique(const int &atom1, const int &atom2, const int &atom3) const;

    void minimiseNetwork();
    double getPotentialEnergy();
    void startMovie();

    void writeMovie();
    void stopMovie();

    void saveState() const;
    void recoverState();
    void setAtomCoords(const int &atomID, const std::vector<double> &newCoords, const int &dim);
};

#endif // LAMMPS_OBJECT_H
