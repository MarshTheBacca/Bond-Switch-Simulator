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

#include "vector_tools.h"

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

struct LammpsObject {
    Network networkA;
    Network networkB;

    void *handle;
    int version;
    int natoms = 0;
    int nbonds = 0;
    int nangles = 0;
    double *bonds = nullptr;

    std::vector<int> angleHelper = std::vector<int>(6);

    LoggerPtr logger;

    LammpsObject();
    explicit LammpsObject(const LoggerPtr &loggerArg);

    void minimiseNetwork();
    double getPotentialEnergy();

    std::vector<double> getCoords(const int &dim) const;
    void setCoords(std::vector<double> &newCoords, int dim);
    void setAtomCoords(const int &atomID, const std::vector<double> &newCoords, const int &dim);

    void breakBond(const int &atom1, const int &atom2, const int &type);
    void formBond(const int &atom1, const int &atom2, const int &type);
    void breakAngle(const int &atom1, const int &atom2, const int &atom3);
    void formAngle(const int &atom1, const int &atom2, const int &atom3);
    void breakAndFormAngles(int atom3, int atom1, int atom5, int atom4, int e1, int e11, int d1, int d11);
    std::pair<int, int> identifyAtoms(int atomA, int atomB, Network networkAArg);

    void switchGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                        const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes,
                        const std::vector<double> &rotatedCoord1, const std::vector<double> &rotatedCoord2);
    void revertGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                        const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes);
    std::vector<int> getAngles() const;


    void writeData();
    void startMovie();
    void writeMovie();
    void stopMovie();

    void showAngles(const int &numLines) const;
    bool checkAngleUnique(const int &atom1, const int &atom2, const int &atom3) const;

};

#endif // LAMMPS_OBJECT_H
