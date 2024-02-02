// Created by olwhi on 24/07/2023, edited by Marshall Hunt 28/01/2024.
#include "lammps_object.h"
#include <stdio.h>
#include <cstdint>
#include <unistd.h>
#include <fstream>

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
LammpsObject::LammpsObject(const std::string &selector, const std::string &inputFolder, const LoggerPtr logger)
{
    logger->info("Creating Lammps Object");
    const char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(const char *);
    logger->info("lmpargc -> {} ", lmpargc);
    handle = lammps_open_no_mpi(lmpargc, const_cast<char **>(lmpargv), nullptr);

    if (handle == nullptr)
    {
        logger->critical("LAMMPS initialization failed");
        lammps_mpi_finalize();
        throw std::runtime_error("LAMMPS initialization failed");
    }
    logger->info("Created LAMMPS handle");
    version = lammps_version(handle);
    logger->info("LAMMPS Version: {}", version);

    std::map<std::string, std::string> selectorToFile = {
        {"Si", "Si.in"},
        {"Si2O3", "Si2O3.in"},
        {"SiO2", "SiO2.in"},
        {"C", "C.in"},
        {"BN", "BN.in"}};

    if (selectorToFile.find(selector) == selectorToFile.end())
    {
        logger->critical("Invalid Selector");
        throw std::runtime_error("Invalid Selector");
    }

    std::string inputFilePath = inputFolder + "/" + selectorToFile[selector];
    logger->info("Executing LAMMPS Script: {}", inputFilePath);
    lammps_file(handle, inputFilePath.c_str());
    natoms = (int)(lammps_get_natoms(handle) + 0.5);
    logger->info("LAMMPS Number of atoms {}", natoms);
    if (const auto nbonds_ptr = static_cast<const int *>(lammps_extract_global(handle, "nbonds")); nbonds_ptr)
    {
        nbonds = *nbonds_ptr;
        logger->info("LAMMPS Number of bonds {}", nbonds);
    }
    else
    {
        logger->error("Failed to extract number of bonds");
    }

    if (const auto nangles_ptr = static_cast<int *>(lammps_extract_global(handle, "nangles")); nangles_ptr)
    {
        nangles = *nangles_ptr;
        logger->info("LAMMPS Number of angles {}", nangles);
    }
    else
    {
        logger->error("Failed to extract number of angles");
    }
}

void LammpsObject::write_data(const std::string &structureName)
{
    std::unordered_map<std::string, std::string> selectorToFile = {
        {"Si", "Si_results.in"},
        {"Si2O3", "Si2O3_results.in"},
        {"SiO2", "SiO2_results.in"},
        {"C", "C_results.in"},
        {"BN", "BN_results.in"}};

    auto iterator = selectorToFile.find(structureName);
    if (iterator == selectorToFile.end())
    {
        throw std::runtime_error("Invalid structure name");
    }
    std::string command = "write_data " + prefixFolderOut + "/" + iterator->second;
    lammps_command(handle, command.c_str());
}

void LammpsObject::write_restart(const std::string &structureName)
{
    std::unordered_map<std::string, std::string> structureToFile = {
        {"Si", "Si_restart.restart"},
        {"Si2O3", "Si2O3_restart.restart"},
        {"SiO2", "SiO2_restart.restart"},
        {"C", "C_restart.restart"},
        {"BN", "BN_restart.restart"}};

    // the iterator is a pointer to a key value pair
    auto iterator = structureToFile.find(structureName);
    if (iterator == structureToFile.end())
    {
        throw std::runtime_error("Invalid structure name");
    }

    // that's why we use iterator->second to access the value of the pair
    std::string command = "write_restart " + prefixFolderOut + "/" + iterator->second;
    lammps_command(handle, command.c_str());
}

void LammpsObject::finaliseLAMMPSObject(const std::string &structureName)
{
    std::unordered_map<std::string, std::string> structureToFile = {
        {"Si", "Si_restart.restart"},
        {"Si2O3", "Si2O3_restart.restart"},
        {"SiO2", "SiO2_restart.restart"},
        {"C", "C_restart.restart"},
        {"BN", "BN_restart.restart"}};
    auto iterator = structureToFile.find(structureName);
    if (iterator == structureToFile.end())
    {
        throw std::runtime_error("Invalid structure name");
    }
    std::string command = "write_data " + prefixFolderOut + "/" + iterator->second;
    lammps_command(handle, command.c_str());
    lammps_close(handle);
}

void LammpsObject::runInput(const std::string &fname)
{
    lammps_file(handle, fname.c_str());
    nbonds = std::intptr_t(lammps_extract_global(handle, "nbonds"));
    nangles = std::intptr_t(lammps_extract_global(handle, "nangles"));
}

double LammpsObject::pbx()
{
    return *(double *)lammps_extract_global(handle, "boxxhi");
}
double LammpsObject::pby()
{
    return *(double *)lammps_extract_global(handle, "boxyhi");
}
double LammpsObject::pbz()
{
    return *(double *)lammps_extract_global(handle, "boxzhi");
}

double *LammpsObject::fetchCrds(int dim)
{
    coords = new double[dim * natoms];
    lammps_gather_atoms(handle, "x", 1, dim, coords);
    return coords;
}

void LammpsObject::pushCrds(int dim, double *old_coords)
{
    lammps_scatter_atoms(handle, "x", 1, dim, old_coords);
}

double *LammpsObject::fetchBonds()
{
    double initialBondCount = lammps_get_thermo(handle, "bonds");
    int bond_size = 3 * initialBondCount;
    std::cout << bond_size << std::endl;
    bonds = new double[bond_size];
    lammps_gather_bonds(handle, bonds);
    return bonds;
}

void LammpsObject::breakBond(int atom1, int atom2, int type, LoggerPtr logger)
{
    logger->debug("Breaking bond between {} and {} of type {}", atom1, atom2, type);
    double initialBondCount = lammps_get_thermo(handle, "bonds");
    std::string command;
    command = "group switch id " + std::to_string(atom1) + " " + std::to_string(atom2);
    lammps_command(handle, command.c_str());
    command = "delete_bonds switch bond " + std::to_string(type) + " remove";
    lammps_command(handle, command.c_str());
    lammps_command(handle, "group switch delete");
    double finalBondCount = lammps_get_thermo(handle, "bonds");
    if (finalBondCount != initialBondCount - 1)
    {
        std::ostringstream oss;
        oss << "Error in Bond Counts while breaking bond, initial: " << finalBondCount << " final: " << initialBondCount;
        logger->critical(oss.str());
        throw std::runtime_error(oss.str());
    }
}

void LammpsObject::formBond(int atom1, int atom2, int type, LoggerPtr logger)
{
    logger->debug("Forming bond between {} and {} of type {}", atom1, atom2, type);
    std::string command;
    command = "create_bonds single/bond " + std::to_string(type) + " " + std::to_string(atom1) + " " + std::to_string(atom2);
    lammps_command(handle, command.c_str());
}

void LammpsObject::breakAngle(int atom1, int atom2, int atom3)
{
    double initialAngleCount = lammps_get_thermo(handle, "angles");

    std::string command;
    command = "group switch id " + std::to_string(atom1) + " " + std::to_string(atom2) + " " + std::to_string(atom3);
    lammps_command(handle, command.c_str());
    lammps_command(handle, "delete_bonds switch angle 1 remove");
    lammps_command(handle, "group switch delete");
    double finalAngleCount = lammps_get_thermo(handle, "angles");
    if (finalAngleCount != initialAngleCount - 1)
    {
        std::cout << "Error in Angle Counts" << std::endl;
        exit(7);
    }
}

void LammpsObject::formAngle(int atom1, int atom2, int atom3)
{
    std::string command;
    command = "create_bonds single/angle 1 " + std::to_string(atom1) + " " + std::to_string(atom2) + " " + std::to_string(atom3);
    lammps_command(handle, command.c_str());
}

// solution here likely involves grabbing angles and bonds from Network
void LammpsObject::switchGraphene(VecF<int> switchIdsA, Network networkA, LoggerPtr logger)
{
    // unpck parameters
    int atom1 = switchIdsA[0] + 1;
    int atom2 = switchIdsA[1] + 1;
    int atom3 = switchIdsA[2] + 1;
    int atom4 = switchIdsA[3] + 1;
    int atom5 = switchIdsA[4] + 1;
    int atom6 = switchIdsA[5] + 1;

    breakBond(atom1, atom5, 1, logger);
    breakBond(atom2, atom4, 1, logger);
    formBond(atom1, atom4, 1, logger);
    formBond(atom2, atom5, 1, logger);

    int e1;
    int e11;
    int d1;
    int d11;
    if (networkA.nodes[atom5 - 1].netCnxs[0] == atom2 - 1)
    {
        e1 = networkA.nodes[atom5 - 1].netCnxs[1] + 1;
        e11 = networkA.nodes[atom5 - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[atom5 - 1].netCnxs[1] == atom2 - 1)
    {
        e1 = networkA.nodes[atom5 - 1].netCnxs[0] + 1;
        e11 = networkA.nodes[atom5 - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[atom5 - 1].netCnxs[2] == atom2 - 1)
    {
        e1 = networkA.nodes[atom5 - 1].netCnxs[0] + 1;
        e11 = networkA.nodes[atom5 - 1].netCnxs[1] + 1;
    }
    else
    {
        logger->warn("atom5 {} not connected to atom2 {}", atom5, atom2);
        sleep(10);
    }
    if (networkA.nodes[atom4 - 1].netCnxs[0] == atom1 - 1)
    {
        d1 = networkA.nodes[atom4 - 1].netCnxs[1] + 1;
        d11 = networkA.nodes[atom4 - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[atom4 - 1].netCnxs[1] == atom1 - 1)
    {
        d1 = networkA.nodes[atom4 - 1].netCnxs[0] + 1;
        d11 = networkA.nodes[atom4 - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[atom4 - 1].netCnxs[2] == atom1 - 1)
    {
        d1 = networkA.nodes[atom4 - 1].netCnxs[0] + 1;
        d11 = networkA.nodes[atom4 - 1].netCnxs[1] + 1;
    }
    else
    {
        logger->warn("atom4 {} not connected to atom1 {}", atom4, atom1);
        sleep(10);
    }

    breakAngle(atom3, atom1, atom5);
    breakAngle(atom2, atom1, atom5);
    breakAngle(atom4, atom2, atom6);
    breakAngle(atom4, atom2, atom1);

    breakAngle(e1, atom5, atom1);
    breakAngle(e11, atom5, atom1);
    breakAngle(d1, atom4, atom2);
    breakAngle(d11, atom4, atom2);

    formAngle(atom3, atom1, atom4);
    formAngle(atom2, atom1, atom4);
    formAngle(atom6, atom2, atom5);
    formAngle(atom1, atom2, atom5);

    formAngle(e1, atom5, atom2);
    formAngle(e11, atom5, atom2);
    formAngle(d1, atom4, atom1);
    formAngle(d11, atom4, atom1);
}

void LammpsObject::revertGraphene(VecF<int> switchIdsA, Network networkA, LoggerPtr logger)
{

    // unpck parameters
    int a = switchIdsA[0] + 1;
    int b = switchIdsA[1] + 1;
    int c = switchIdsA[2] + 1;
    int d = switchIdsA[3] + 1;
    int e = switchIdsA[4] + 1;
    int f = switchIdsA[5] + 1;

    breakBond(a, d, 1, logger);
    breakBond(b, e, 1, logger);
    formBond(a, e, 1, logger);
    formBond(b, d, 1, logger);

    int e1;
    int e11;
    int d1;
    int d11;
    if (networkA.nodes[e - 1].netCnxs[0] == a - 1)
    {
        e1 = networkA.nodes[e - 1].netCnxs[1] + 1;
        e11 = networkA.nodes[e - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[e - 1].netCnxs[1] == a - 1)
    {
        e1 = networkA.nodes[e - 1].netCnxs[0] + 1;
        e11 = networkA.nodes[e - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[e - 1].netCnxs[2] == a - 1)
    {
        e1 = networkA.nodes[e - 1].netCnxs[0] + 1;
        e11 = networkA.nodes[e - 1].netCnxs[1] + 1;
    }
    else
    {
        logger->warn("e {} not connected to a {}, sleeping 10", e, a);
        sleep(10);
    }
    if (networkA.nodes[d - 1].netCnxs[0] == b - 1)
    {
        d1 = networkA.nodes[d - 1].netCnxs[1] + 1;
        d11 = networkA.nodes[d - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[d - 1].netCnxs[1] == b - 1)
    {
        d1 = networkA.nodes[d - 1].netCnxs[0] + 1;
        d11 = networkA.nodes[d - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[d - 1].netCnxs[2] == b - 1)
    {
        d1 = networkA.nodes[d - 1].netCnxs[0] + 1;
        d11 = networkA.nodes[d - 1].netCnxs[1] + 1;
    }
    else
    {
        logger->warn("d {} not connected to b {}, sleeping 10", d, b);
        sleep(10);
    }

    breakAngle(c, a, d);
    breakAngle(b, a, d);
    breakAngle(f, b, e);
    breakAngle(a, b, e);

    breakAngle(e1, e, b);
    breakAngle(e11, e, b);
    breakAngle(d1, d, a);
    breakAngle(d11, d, a);

    formAngle(c, a, e);
    formAngle(b, a, e);
    formAngle(d, b, f);
    formAngle(d, b, a);

    formAngle(e1, e, a);
    formAngle(e11, e, a);
    formAngle(d1, d, b);
    formAngle(d11, d, b);
}

void LammpsObject::switchTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger)
{
    // unpck parameters
    int a = switchIdsA[0] + 1;
    int b = switchIdsA[1] + 1;

    int beta = switchIdsT[7] + 1;
    int gamma = switchIdsT[8] + 1;
    int delta = switchIdsT[9] + 1;
    int eta = switchIdsT[10] + 1;

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

void LammpsObject::revertTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger)
{
    // unpck parameters
    int a = switchIdsA[0] + 1;
    int b = switchIdsA[1] + 1;

    int beta = switchIdsT[7] + 1;
    int gamma = switchIdsT[8] + 1;
    int delta = switchIdsT[9] + 1;
    int eta = switchIdsT[10] + 1;

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

void LammpsObject::switchBilayer(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger)
{

    // unpck parameters
    int a = switchIdsA[0];
    int b = switchIdsA[1];

    int beta = switchIdsT[7];
    int gamma = switchIdsT[8];
    int delta = switchIdsT[9];
    int eta = switchIdsT[10];

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

void LammpsObject::revertBilayer(VecF<int> switchIdsA, VecF<int> switchIdsT, LoggerPtr logger)
{

    // unpck parameters
    int a = switchIdsA[0];
    int b = switchIdsA[1];

    int beta = switchIdsT[7];
    int gamma = switchIdsT[8];
    int delta = switchIdsT[9];
    int eta = switchIdsT[10];

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
    int nO;

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

void LammpsObject::switchBonds(VecF<int> switchIdsA, VecF<int> switchIdsT)
{
    // unpck parameters
    int a = switchIdsA[0];
    int b = switchIdsA[1];
    int d = switchIdsA[3];
    int e = switchIdsA[4];

    int beta = switchIdsT[7];
    int gamma = switchIdsT[8];
    int delta = switchIdsT[9];
    int eta = switchIdsT[10];

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

    int nSi;
    nSi = round(natoms / 3);

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

    a_prime = 2 * a;
    a_prime_prime = 2 * a + 1;
    a_prime_prime_prime = 3 * a;

    b_prime = 2 * b;
    b_prime_prime = 2 * b + 1;
    b_prime_prime_prime = 3 * b;

    beta_prime = 3 * nSi + (beta - nSi) * 2;
    beta_prime_prime = beta_prime + 1;
    gamma_prime = 3 * nSi + (gamma - nSi) * 2;
    gamma_prime_prime = gamma_prime + 1;
    eta_prime = 3 * nSi + (eta - nSi) * 2;
    eta_prime_prime = eta_prime + 1;
    delta_prime = 3 * nSi + (delta - nSi) * 2;
    delta_prime_prime = delta_prime + 1;

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

void LammpsObject::revertBonds(VecF<int> switchIdsA, VecF<int> switchIdsT)
{
    // unpck parameters
    int a = switchIdsA[0];
    int b = switchIdsA[1];
    int d = switchIdsA[3];
    int e = switchIdsA[4];

    int beta = switchIdsT[7];
    int gamma = switchIdsT[8];
    int delta = switchIdsT[9];
    int eta = switchIdsT[10];
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
    command = "group switch id " + std::to_string(a) + " " + std::to_string(d);
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

VecF<int> LammpsObject::GlobalPotentialMinimisation()
{
    lammps_command(handle, "minimize 1.0e-6 0.0 1000000 10000000");
    VecF<int> optstatus(2);
    optstatus[0] = 0;
    optstatus[1] = 10;
    return optstatus;
}

double LammpsObject::GlobalPotentialEnergy()
{
    return lammps_get_thermo(handle, "pe");
}

double LammpsObject::probeE()
{
    return lammps_get_thermo(handle, "etotal");
}