// Created by olwhi on 24/07/2023, edited by Marshall Hunt 28/01/2024.
#include "lammps_object.h"
#include <cstdint>
#include <fstream>
#include <stdio.h>
#include <unistd.h>

/**
 * @brief Default constructor for a blank Lammps Object
 */
LammpsObject::LammpsObject() = default;

/**
 * @brief Constructor for a Lammps Object
 * @param selector The selector for the structure to be created
 * @param inputFolder The folder containing the input files
 * @param logger The logger object
 */
LammpsObject::LammpsObject(const std::string &structureName, const std::string &inputFolder, const LoggerPtr logger) {
    logger->info("Creating Lammps Object");
    const char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(const char *);
    logger->info("lmpargc -> {} ", lmpargc);
    handle = lammps_open_no_mpi(lmpargc, const_cast<char **>(lmpargv), nullptr);

    if (handle == nullptr) {
        logger->critical("LAMMPS initialization failed");
        lammps_mpi_finalize();
        throw std::runtime_error("LAMMPS initialization failed");
    }
    logger->debug("Created LAMMPS handle");
    version = lammps_version(handle);
    logger->debug("LAMMPS Version: {}", version);

    std::map<std::string, std::string> selectorToFile = {
        {"Si", "Si.in"},
        {"Si2O3", "Si2O3.in"},
        {"SiO2", "SiO2.in"},
        {"C", "C.in"},
        {"BN", "BN.in"}};

    if (selectorToFile.find(structureName) == selectorToFile.end()) {
        logger->critical("Invalid Selector");
        throw std::runtime_error("Invalid Selector");
    }

    std::string inputFilePath = inputFolder + "/" + selectorToFile[structureName];
    logger->info("Executing LAMMPS Script: {}", inputFilePath);
    lammps_file(handle, inputFilePath.c_str());
    natoms = (int)(lammps_get_natoms(handle) + 0.5);
    if (const auto nbonds_ptr = static_cast<const int *>(lammps_extract_global(handle, "nbonds")); nbonds_ptr) {
        nbonds = *nbonds_ptr;
    } else {
        logger->error("Failed to extract number of bonds");
    }

    if (const auto nangles_ptr = static_cast<int *>(lammps_extract_global(handle, "nangles")); nangles_ptr) {
        nangles = *nangles_ptr;
    } else {
        logger->error("Failed to extract number of angles");
    }
    logger->info("LAMMPS #nodes: {} #bonds: {} #angles: {}", natoms, nbonds, nangles);
}

void LammpsObject::write_data(const std::string &structureName) {
    std::unordered_map<std::string, std::string> selectorToFile = {
        {"Si", "Si_results.in"},
        {"Si2O3", "Si2O3_results.in"},
        {"SiO2", "SiO2_results.in"},
        {"C", "C_results.in"},
        {"BN", "BN_results.in"}};

    auto iterator = selectorToFile.find(structureName);
    if (iterator == selectorToFile.end()) {
        throw std::runtime_error("Invalid structure name");
    }
    std::string command = "write_data " + prefixFolderOut + "/" + iterator->second;
    lammps_command(handle, command.c_str());
}

void LammpsObject::write_restart(const std::string &structureName) {
    std::unordered_map<std::string, std::string> structureToFile = {
        {"Si", "Si_restart.restart"},
        {"Si2O3", "Si2O3_restart.restart"},
        {"SiO2", "SiO2_restart.restart"},
        {"C", "C_restart.restart"},
        {"BN", "BN_restart.restart"}};

    // the iterator is a pointer to a key value pair
    auto iterator = structureToFile.find(structureName);
    if (iterator == structureToFile.end()) {
        throw std::runtime_error("Invalid structure name");
    }

    // that's why we use iterator->second to access the value of the pair
    std::string command = "write_restart " + prefixFolderOut + "/" + iterator->second;
    lammps_command(handle, command.c_str());
}

void LammpsObject::finaliseLAMMPSObject(const std::string &structureName) {
    std::unordered_map<std::string, std::string> structureToFile = {
        {"Si", "Si_restart.restart"},
        {"Si2O3", "Si2O3_restart.restart"},
        {"SiO2", "SiO2_restart.restart"},
        {"C", "C_restart.restart"},
        {"BN", "BN_restart.restart"}};
    auto iterator = structureToFile.find(structureName);
    if (iterator == structureToFile.end()) {
        throw std::runtime_error("Invalid structure name");
    }
    std::string command = "write_data " + prefixFolderOut + "/" + iterator->second;
    lammps_command(handle, command.c_str());
    lammps_close(handle);
}

double LammpsObject::pbx() {
    return *(double *)lammps_extract_global(handle, "boxxhi");
}
double LammpsObject::pby() {
    return *(double *)lammps_extract_global(handle, "boxyhi");
}
double LammpsObject::pbz() {
    return *(double *)lammps_extract_global(handle, "boxzhi");
}

double *LammpsObject::fetchCrds(int dim) {
    coords = new double[dim * natoms];
    lammps_gather_atoms(handle, "x", 1, dim, coords);
    return coords;
}

void LammpsObject::pushCrds(int dim, double *old_coords) {
    lammps_scatter_atoms(handle, "x", 1, dim, old_coords);
}

void LammpsObject::breakBond(int atom1, int atom2, int type, LoggerPtr logger) {
    logger->debug("Breaking bond between {} and {} of type {}", atom1, atom2, type);
    double initialBondCount = lammps_get_thermo(handle, "bonds");

    // Compute bond properties
    std::string command = "compute myBonds all property/local batom1 batom2 btype";
    lammps_command(handle, command.c_str());

    // Output bond properties
    command = "dump myDump all local 1 bond.dump index c_myBonds[1] c_myBonds[2] c_myBonds[3]";
    lammps_command(handle, command.c_str());

    command = "group switch id " + std::to_string(atom1) + " " + std::to_string(atom2);
    lammps_command(handle, command.c_str());

    // Delete the bond
    command = "delete_bonds switch bond " + std::to_string(type) + " remove";
    lammps_command(handle, command.c_str());

    // Stop outputting bond properties
    command = "undump myDump";
    lammps_command(handle, command.c_str());

    lammps_command(handle, "group switch delete");

    command = "uncompute myBonds";
    lammps_command(handle, command.c_str());

    double finalBondCount = lammps_get_thermo(handle, "bonds");
    if (finalBondCount != initialBondCount - 1) {
        std::ostringstream oss;
        oss << "Error in Bond Counts while breaking bond between " << atom1 << " and " << atom2 << ", initial: " << finalBondCount << " final: " << initialBondCount;
        throw std::runtime_error(oss.str());
    }
}

void LammpsObject::formBond(int atom1, int atom2, int type, LoggerPtr logger) {
    logger->debug("Forming bond between {} and {} of type {}", atom1, atom2, type);
    double initialBondCount = lammps_get_thermo(handle, "bonds");
    std::string command;
    command = "create_bonds single/bond " + std::to_string(type) + " " + std::to_string(atom1) + " " + std::to_string(atom2);
    lammps_command(handle, command.c_str());
    double finalBondCount = lammps_get_thermo(handle, "bonds");
    if (finalBondCount != initialBondCount + 1) {
        std::ostringstream oss;
        oss << "Error in Bond Counts while forming bond, initial: " << finalBondCount << " final: " << initialBondCount;
        throw std::runtime_error(oss.str());
    }
}

/**
 * @brief Breaks an angle in the lattice by trying 1-2-3 first, then 3-2-1
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 * @param logger The logger object
 * @throws std::runtime_error if the angle count doesn't decrease
 */
void LammpsObject::breakAngle(int atom1, int atom2, int atom3, LoggerPtr logger) {
    logger->debug("Breaking angle between {}, {}, {}", atom1, atom2, atom3);
    double initialAngleCount = lammps_get_thermo(handle, "angles");

    // Try to break the angle in both directions
    for (int attempt = 0; attempt < 2; ++attempt) {
        std::string command;
        // Define a group, "switch", containing the atoms in the angle
        if (attempt == 0)
            command = "group switch id " + std::to_string(atom1) + " " + std::to_string(atom2) + " " + std::to_string(atom3);
        else
            command = "group switch id " + std::to_string(atom3) + " " + std::to_string(atom2) + " " + std::to_string(atom1);
        lammps_command(handle, command.c_str());
        lammps_command(handle, "delete_bonds switch angle 1 remove");

        // Delete the user defined group
        lammps_command(handle, "group switch delete");

        // Check if the number of angles has decreased
        double finalAngleCount = lammps_get_thermo(handle, "angles");
        if (finalAngleCount == initialAngleCount - 1)
            return; // The angle was successfully broken
    }

    // If the function hasn't returned by now, the angle couldn't be broken in either direction
    std::stringstream error;
    error << "Error in Angle Counts, initial: " << initialAngleCount;
    logger->critical(error.str());
    throw std::runtime_error(error.str());
}

/**
 * @brief Forms an angle in the lattice
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 * @param logger The logger object
 */
void LammpsObject::formAngle(int atom1, int atom2, int atom3, LoggerPtr logger) {
    logger->debug("Forming angle between {}, {}, {}", atom1, atom2, atom3);
    double initialAngleCount = lammps_get_thermo(handle, "angles");
    std::string command;
    command = "create_bonds single/angle 1 " + std::to_string(atom1) + " " + std::to_string(atom2) + " " + std::to_string(atom3);
    lammps_command(handle, command.c_str());
    double finalAngleCount = lammps_get_thermo(handle, "angles");
    if (finalAngleCount != initialAngleCount + 1) {
        std::ostringstream oss;
        oss << "Error in Angle Counts, initial: " << finalAngleCount << " final: " << initialAngleCount;
        throw std::runtime_error(oss.str());
    }
}

/**
 * @brief perform bond switch in the lattice using NetMC node IDs, ie, zero indexed
 * @param switchIDsA The node IDs of the atoms to be switched
 * @param logger The logger object
 */
void LammpsObject::switchGraphene(VecF<int> switchIDsA, LoggerPtr logger) {
    /*
     *                7-----8                               7-----8
     *               /       \                              |     |
     *              /         \                      11-----3  2  4-----12
     *      11-----3     2     4-----12                      \   /
     *              \         /                               \ /
     *               \       /                                 1
     *          3     1-----2     4         --->         3     |      4
     *               /       \                                 2
     *              /         \                               /  \
     *      13-----5     1     6-----14                      /    \
     *              \         /                      13-----5  1   6-----14
     *               \       /                              |      |
     *                9-----10                              9------10
     *
     *      Bonds to break       Bonds to Make
     *      1-5, 2-4             1-4, 2-5
     *
     *      Angles to break      Angles to Make
     *      1-5-9, 1-5-13        1-4-8, 1-4-12
     *      2-4-8, 2-4-12        2-5-9, 2-5-13
     *      4-2-1, 4-2-6         4-1-2, 4-1-3
     *      5-1-2, 5-1-3         6-2-1, 6-2-5
     */

    int atom1 = switchIDsA[0] + 1;
    int atom2 = switchIDsA[1] + 1;
    int atom5 = switchIDsA[2] + 1;
    int atom6 = switchIDsA[3] + 1;
    int atom3 = switchIDsA[4] + 1;
    int atom4 = switchIDsA[5] + 1;
    int atom7 = switchIDsA[6] + 1;
    int atom8 = switchIDsA[7] + 1;
    int atom9 = switchIDsA[8] + 1;
    int atom10 = switchIDsA[9] + 1;
    int atom11 = switchIDsA[10] + 1;
    int atom12 = switchIDsA[11] + 1;
    int atom13 = switchIDsA[12] + 1;
    int atom14 = switchIDsA[13] + 1;
    logger->debug("");
    logger->debug("           {}----{}                                {}------{}", atom7, atom8, atom7, atom8);
    logger->debug("          /        \\                                |      |");
    logger->debug("         /          \\                       {}-----{}     {}-----{}", atom11, atom3, atom4, atom12);
    logger->debug(" {}-----{}          {}-----{}                        \\    /", atom11, atom3, atom4, atom12);
    logger->debug("         \\          /                                 \\  /");
    logger->debug("          \\        /                                   {}", atom1);
    logger->debug("          {}-----{}                --->                |       ", atom1, atom2);
    logger->debug("          /       \\                                    {}", atom2);
    logger->debug("         /         \\                                  /  \\");
    logger->debug(" {}-----{}         {}-----{}                         /    \\", atom13, atom5, atom6, atom14);
    logger->debug("         \\         /                        {}-----{}     {}-----{}", atom13, atom5, atom6, atom14);
    logger->debug("          \\       /                                 |      |");
    logger->debug("           {}----{}                                {}------{}", atom9, atom10, atom9, atom10);
    logger->debug("");

    breakBond(atom1, atom5, 1, logger);
    breakBond(atom2, atom4, 1, logger);
    formBond(atom1, atom4, 1, logger);
    formBond(atom2, atom5, 1, logger);

    breakAngle(atom1, atom5, atom9, logger);
    breakAngle(atom1, atom5, atom13, logger);
    breakAngle(atom2, atom4, atom8, logger);
    breakAngle(atom2, atom4, atom12, logger);
    breakAngle(atom4, atom2, atom1, logger);
    breakAngle(atom4, atom2, atom6, logger);
    breakAngle(atom5, atom1, atom2, logger);
    breakAngle(atom5, atom1, atom3, logger);

    formAngle(atom1, atom4, atom8, logger);
    formAngle(atom1, atom4, atom12, logger);
    formAngle(atom2, atom5, atom9, logger);
    formAngle(atom2, atom5, atom13, logger);
    formAngle(atom4, atom1, atom2, logger);
    formAngle(atom4, atom1, atom3, logger);
    formAngle(atom5, atom2, atom1, logger);
    formAngle(atom6, atom2, atom5, logger);
}

void LammpsObject::revertGraphene(VecF<int> switchIDsA, LoggerPtr logger) {
    int atom1 = switchIDsA[0] + 1;
    int atom2 = switchIDsA[1] + 1;
    int atom5 = switchIDsA[2] + 1;
    int atom6 = switchIDsA[3] + 1;
    int atom3 = switchIDsA[4] + 1;
    int atom4 = switchIDsA[5] + 1;
    int atom7 = switchIDsA[6] + 1;
    int atom8 = switchIDsA[7] + 1;
    int atom9 = switchIDsA[8] + 1;
    int atom10 = switchIDsA[9] + 1;
    int atom11 = switchIDsA[10] + 1;
    int atom12 = switchIDsA[11] + 1;
    int atom13 = switchIDsA[12] + 1;
    int atom14 = switchIDsA[13] + 1;

    breakBond(atom1, atom4, 1, logger);
    breakBond(atom2, atom5, 1, logger);
    formBond(atom1, atom5, 1, logger);
    formBond(atom2, atom4, 1, logger);

    breakAngle(atom1, atom4, atom8, logger);
    breakAngle(atom1, atom4, atom12, logger);
    breakAngle(atom2, atom5, atom9, logger);
    breakAngle(atom2, atom5, atom13, logger);
    breakAngle(atom4, atom1, atom2, logger);
    breakAngle(atom4, atom1, atom3, logger);
    breakAngle(atom5, atom2, atom1, logger);
    breakAngle(atom6, atom2, atom5, logger);

    formAngle(atom1, atom5, atom9, logger);
    formAngle(atom1, atom5, atom13, logger);
    formAngle(atom2, atom4, atom8, logger);
    formAngle(atom2, atom4, atom12, logger);
    formAngle(atom4, atom2, atom1, logger);
    formAngle(atom4, atom2, atom6, logger);
    formAngle(atom5, atom1, atom2, logger);
    formAngle(atom5, atom1, atom3, logger);
}

void LammpsObject::switchTriangleRaft(VecF<int> switchIDsA, VecF<int> switchIDsT, LoggerPtr logger) {
    // unpck parameters
    int a = switchIDsA[0] + 1;
    int b = switchIDsA[1] + 1;

    int beta = switchIDsT[7] + 1;
    int gamma = switchIDsT[8] + 1;
    int delta = switchIDsT[9] + 1;
    int eta = switchIDsT[10] + 1;

    // Corrections to
    //  -- Actually remove the bonds
    //  -- Correct bond types

    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */
    breakBond(a, gamma, 2, logger);
    breakBond(b, delta, 2, logger);
    breakBond(gamma, beta, 1, logger);
    breakBond(delta, eta, 1, logger);
    formBond(a, delta, 2, logger);
    formBond(b, gamma, 2, logger);
    formBond(beta, delta, 1, logger);
    formBond(eta, gamma, 1, logger);
}

void LammpsObject::revertTriangleRaft(VecF<int> switchIDsA, VecF<int> switchIDsT, LoggerPtr logger) {
    // unpck parameters
    int a = switchIDsA[0] + 1;
    int b = switchIDsA[1] + 1;

    int beta = switchIDsT[7] + 1;
    int gamma = switchIDsT[8] + 1;
    int delta = switchIDsT[9] + 1;
    int eta = switchIDsT[10] + 1;

    breakBond(a, delta, 2, logger);
    breakBond(b, gamma, 2, logger);
    breakBond(beta, delta, 1, logger);
    breakBond(eta, gamma, 1, logger);
    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */

    formBond(a, gamma, 2, logger);
    formBond(b, delta, 2, logger);
    formBond(gamma, beta, 1, logger);
    formBond(delta, eta, 1, logger);
}

void LammpsObject::switchBilayer(VecF<int> switchIDsA, VecF<int> switchIDsT, LoggerPtr logger) {

    // unpck parameters
    int a = switchIDsA[0];
    int b = switchIDsA[1];

    int beta = switchIDsT[7];
    int gamma = switchIDsT[8];
    int delta = switchIDsT[9];
    int eta = switchIDsT[10];

    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */

    int nSi;
    int nSi2;
    nSi = (int)(round(natoms / 3) + 0.5);
    nSi2 = (int)(round(nSi / 2) + 0.5);

    int a_prime;
    int a_prime_prime;
    int b_prime;
    int b_prime_prime;
    int a_prime_prime_prime;
    int b_prime_prime_prime;

    int beta_prime;
    int beta_prime_prime;
    int gamma_prime;
    int gamma_prime_prime;
    int eta_prime;
    int eta_prime_prime;
    int delta_prime;
    int delta_prime_prime;

    a_prime = 2 * a + 1;
    a_prime_prime = 2 * a + 2;
    a_prime_prime_prime = nSi + a + 1;

    b_prime = 2 * b + 1;
    b_prime_prime = 2 * b + 2;
    b_prime_prime_prime = nSi + b + 1;

    beta_prime = 3 * nSi2 + (beta - nSi2) * 2 + 1;
    beta_prime_prime = beta_prime + 1;
    gamma_prime = 3 * nSi2 + (gamma - nSi2) * 2 + 1;
    gamma_prime_prime = gamma_prime + 1;
    eta_prime = 3 * nSi2 + (eta - nSi2) * 2 + 1;
    eta_prime_prime = eta_prime + 1;
    delta_prime = 3 * nSi2 + (delta - nSi2) * 2 + 1;
    delta_prime_prime = delta_prime + 1;

    //      Top layer
    breakBond(a_prime, gamma_prime, 2, logger);
    breakBond(b_prime, delta_prime, 2, logger);
    breakBond(beta_prime, gamma_prime, 1, logger);
    breakBond(eta_prime, delta_prime, 1, logger);
    formBond(a_prime, delta_prime, 2, logger);
    formBond(b_prime, gamma_prime, 2, logger);
    formBond(delta_prime, beta_prime, 1, logger);
    formBond(gamma_prime, eta_prime, 1, logger);
    // Bottom layer
    breakBond(a_prime_prime, gamma_prime_prime, 2, logger);
    breakBond(b_prime_prime, delta_prime_prime, 2, logger);
    breakBond(beta_prime_prime, gamma_prime_prime, 1, logger);
    breakBond(eta_prime_prime, delta_prime_prime, 1, logger);
    formBond(a_prime_prime, delta_prime_prime, 2, logger);
    formBond(b_prime_prime, gamma_prime_prime, 2, logger);
    formBond(beta_prime_prime, delta_prime_prime, 1, logger);
    formBond(eta_prime_prime, gamma_prime_prime, 1, logger);
    // Bridge to Bottom Layer
    breakBond(a_prime_prime_prime, gamma_prime_prime, 1, logger);
    breakBond(b_prime_prime_prime, delta_prime_prime, 1, logger);

    formBond(a_prime_prime_prime, delta_prime_prime, 1, logger);
    formBond(b_prime_prime_prime, gamma_prime_prime, 1, logger);
    breakBond(a_prime_prime_prime, gamma_prime, 1, logger);
    breakBond(b_prime_prime_prime, delta_prime, 1, logger);
    formBond(a_prime_prime_prime, delta_prime, 1, logger);
    formBond(b_prime_prime_prime, gamma_prime, 1, logger);
}

void LammpsObject::revertBilayer(VecF<int> switchIDsA, VecF<int> switchIDsT, LoggerPtr logger) {

    // unpck parameters
    int a = switchIDsA[0];
    int b = switchIDsA[1];

    int beta = switchIDsT[7];
    int gamma = switchIDsT[8];
    int delta = switchIDsT[9];
    int eta = switchIDsT[10];

    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */

    auto nSi = static_cast<int>(round(natoms / 3) + 0.5);
    auto nSi2 = static_cast<int>(round(nSi / 2) + 0.5);

    int a_prime = 2 * a + 1;
    int a_prime_prime = 2 * a + 2;
    int a_prime_prime_prime = nSi + a + 1;

    int b_prime = 2 * b + 1;
    int b_prime_prime = 2 * b + 2;
    int b_prime_prime_prime = nSi + b + 1;

    int beta_prime = 3 * nSi2 + (beta - nSi2) * 2 + 1;
    int beta_prime_prime = beta_prime + 1;
    int gamma_prime = 3 * nSi2 + (gamma - nSi2) * 2 + 1;
    int gamma_prime_prime = gamma_prime + 1;
    int eta_prime = 3 * nSi2 + (eta - nSi2) * 2 + 1;
    int eta_prime_prime = eta_prime + 1;
    int delta_prime = 3 * nSi2 + (delta - nSi2) * 2 + 1;
    int delta_prime_prime = delta_prime + 1;

    breakBond(a_prime, delta_prime, 2, logger);
    breakBond(b_prime, gamma_prime, 2, logger);
    breakBond(delta_prime, beta_prime, 1, logger);
    breakBond(gamma_prime, eta_prime, 1, logger);
    formBond(a_prime, gamma_prime, 2, logger);
    formBond(b_prime, delta_prime, 2, logger);
    formBond(beta_prime, gamma_prime, 1, logger);
    formBond(eta_prime, delta_prime, 1, logger);
    breakBond(a_prime_prime, delta_prime_prime, 2, logger);
    breakBond(b_prime_prime, gamma_prime_prime, 2, logger);
    breakBond(delta_prime_prime, beta_prime_prime, 1, logger);
    breakBond(gamma_prime_prime, eta_prime_prime, 1, logger);
    formBond(a_prime_prime, gamma_prime_prime, 2, logger);
    formBond(b_prime_prime, delta_prime_prime, 2, logger);
    formBond(beta_prime_prime, gamma_prime_prime, 1, logger);
    formBond(eta_prime_prime, delta_prime_prime, 1, logger);
    breakBond(a_prime_prime_prime, delta_prime, 1, logger);
    breakBond(b_prime_prime_prime, gamma_prime, 1, logger);
    formBond(a_prime_prime_prime, gamma_prime, 1, logger);
    formBond(b_prime_prime_prime, delta_prime, 1, logger);
    breakBond(a_prime_prime_prime, delta_prime_prime, 1, logger);
    breakBond(b_prime_prime_prime, gamma_prime_prime, 1, logger);
    formBond(a_prime_prime_prime, gamma_prime_prime, 1, logger);
    formBond(b_prime_prime_prime, delta_prime_prime, 1, logger);
}

void LammpsObject::switchBonds(VecF<int> switchIDsA, VecF<int> switchIDsT) {
    // unpck parameters
    int a = switchIDsA[0];
    int b = switchIDsA[1];
    int d = switchIDsA[3];
    int e = switchIDsA[4];

    int beta = switchIDsT[7];
    int gamma = switchIDsT[8];
    int delta = switchIDsT[9];
    int eta = switchIDsT[10];

    /* Switch connectivities in lattice and dual
     * 3-3 coordination connection
     * a,b,c,d,e,f are nodes in lattice A
     * u,v,w,x are nodes in lattice B
     *  E      F            V
     *   \    /           / | \
     *   A---B           W  |  X
     *  /     \           \ | /
     * C       D            U
     *
     *       A x E       A - D
     *       B x D       B - E
     *       D x B       D - A
     *       E x A       E - B
     *
     *       Break A - E and B - D
     *       Form  A - D and B - E
     *
     */

    std::string command;
    command = "group switch id " + std::to_string(a) + " " + std::to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b) + " " + std::to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a) + " " + std::to_string(d);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b) + " " + std::to_string(e);
    lammps_command(handle, command.c_str());

    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */
    command = "group switch id " + std::to_string(a) + " " + std::to_string(gamma);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b) + " " + std::to_string(delta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(gamma) + " " + std::to_string(beta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(delta) + " " + std::to_string(eta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a) + " " + std::to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b) + " " + std::to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(beta) + " " + std::to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(eta) + " " + std::to_string(gamma);
    lammps_command(handle, command.c_str());

    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */

    auto nSi = static_cast<int>(round(natoms / 3) + 0.5);

    int a_prime = 2 * a;
    int a_prime_prime = 2 * a + 1;
    int a_prime_prime_prime = 3 * a;

    int b_prime = 2 * b;
    int b_prime_prime = 2 * b + 1;
    int b_prime_prime_prime = 3 * b;

    int beta_prime = 3 * nSi + (beta - nSi) * 2;
    int beta_prime_prime = beta_prime + 1;
    int gamma_prime = 3 * nSi + (gamma - nSi) * 2;
    int gamma_prime_prime = gamma_prime + 1;
    int eta_prime = 3 * nSi + (eta - nSi) * 2;
    int eta_prime_prime = eta_prime + 1;
    int delta_prime = 3 * nSi + (delta - nSi) * 2;
    int delta_prime_prime = delta_prime + 1;

    //      Top layer

    command = "group switch id " + std::to_string(a_prime) + " " + std::to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b_prime) + " " + std::to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(gamma_prime) + " " + std::to_string(beta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(delta_prime) + " " + std::to_string(eta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a_prime) + " " + std::to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b_prime) + " " + std::to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(beta_prime) + " " + std::to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(eta_prime) + " " + std::to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    // Bottom layer

    command = "group switch id " + std::to_string(a_prime_prime) + " " + std::to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b_prime_prime) + " " + std::to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(gamma_prime_prime) + " " + std::to_string(beta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(delta_prime_prime) + " " + std::to_string(eta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a_prime_prime) + " " + std::to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b_prime_prime) + " " + std::to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(beta_prime_prime) + " " + std::to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(eta_prime_prime) + " " + std::to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());

    // Bridge to Top Layer

    command = "group switch id " + std::to_string(a_prime_prime_prime) + " " + std::to_string(gamma_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b_prime_prime_prime) + " " + std::to_string(delta_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a_prime_prime_prime) + " " + std::to_string(delta_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b_prime_prime_prime) + " " + std::to_string(gamma_prime);
    lammps_command(handle, command.c_str());

    // Bridge to Bottom Layer

    command = "group switch id " + std::to_string(a_prime_prime_prime) + " " + std::to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b_prime_prime_prime) + " " + std::to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a_prime_prime_prime) + " " + std::to_string(delta_prime_prime);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b_prime_prime_prime) + " " + std::to_string(gamma_prime_prime);
    lammps_command(handle, command.c_str());
}

void LammpsObject::revertBonds(VecF<int> switchIDsA, VecF<int> switchIDsT) {
    // unpck parameters
    int a = switchIDsA[0];
    int b = switchIDsA[1];
    int d = switchIDsA[3];
    int e = switchIDsA[4];

    int beta = switchIDsT[7];
    int gamma = switchIDsT[8];
    int delta = switchIDsT[9];
    int eta = switchIDsT[10];
    /* Switch connectivities in lattice and dual
     * 3-3 coordination connection
     * a,b,c,d,e,f are nodes in lattice A
     * u,v,w,x are nodes in lattice B
     *  E      F            V
     *   \    /           / | \
     *   A---B           W  |  X
     *  /     \           \ | /
     * C       D            U
     *
     *       A x E       A - D
     *       B x D       B - E
     *       D x B       D - A
     *       E x A       E - B
     *
     *       Break A - E and B - D
     *       Form  A - D and B - E
     *
     */

    std::string command = "group switch id " + std::to_string(a) + " " + std::to_string(d);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b) + " " + std::to_string(e);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a) + " " + std::to_string(e);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b) + " " + std::to_string(d);
    lammps_command(handle, command.c_str());

    /*
     *
     *      A x gamma       A - delta
     *      B x delta       B - gamma
     *
     *      gamma x beta    gamma - eta
     *      beta  x gammma  beta - delta
     *      delta x eta     delta - beta
     *      eta   x delta   eta - gamma
     *
     */
    command = "group switch id " + std::to_string(a) + " " + std::to_string(delta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(b) + " " + std::to_string(gamma);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(gamma) + " " + std::to_string(eta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "group switch id " + std::to_string(delta) + " " + std::to_string(beta);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");

    command = "create_bonds single/bond 1 " + std::to_string(a) + " " + std::to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(b) + " " + std::to_string(delta);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(beta) + " " + std::to_string(gamma);
    lammps_command(handle, command.c_str());

    command = "create_bonds single/bond 1 " + std::to_string(eta) + " " + std::to_string(delta);
    lammps_command(handle, command.c_str());
}

VecF<int> LammpsObject::globalPotentialMinimisation() {
    lammps_command(handle, "minimize 1.0e-6 0.0 1000000 10000000");
    VecF<int> optstatus(2);
    optstatus[0] = 0;
    optstatus[1] = 10;
    return optstatus;
}

double LammpsObject::globalPotentialEnergy() {
    return lammps_get_thermo(handle, "pe");
}

double LammpsObject::probeE() {
    return lammps_get_thermo(handle, "etotal");
}