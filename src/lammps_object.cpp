// Created by olwhi on 24/07/2023, edited by Marshall Hunt 28/01/2024.
#include "lammps_object.h"
#include <stdio.h>
#include <cstdint>
#include <unistd.h>
#include <fstream>

LammpsObject::LammpsObject()
{
}

LammpsObject::LammpsObject(std::string selector, std::string inputFolder, LoggerPtr logger)
{
    logger->info("Creating Lammps Object");
    const char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(const char *);
    logger->info("lmpargc -> {} ", lmpargc);
    handle = lammps_open_no_mpi(lmpargc, (char **)lmpargv, NULL);

    if (handle == NULL)
    {
        logger->critical("LAMMPS initialization failed");
        lammps_mpi_finalize();
        throw std::runtime_error("LAMMPS initialization failed");
    }
    logger->info("Created LAMMPS handle");
    version = lammps_version(handle);
    logger->info("LAMMPS Version: {}", version);

    std::map<std::string, std::string> selectorToFile = {
        {"SimpleGraphene", "/Si.in"},
        {"TriangleRaft", "/Si2O3.in"},
        {"Bilayer", "/SiO2.in"},
        {"Tersoff", "/C.in"},
        {"BN", "/BN.in"}};

    if (selectorToFile.find(selector) == selectorToFile.end())
    {
        logger->critical("Invalid Selector");
        throw std::runtime_error("Invalid Selector");
    }

    std::string inputFilePath = inputFolder + selectorToFile[selector];
    logger->info("Selected {}", selector);
    logger->info("File to run from : {}", inputFilePath);
    lammps_file(handle, inputFilePath.c_str());
    natoms = (int)(lammps_get_natoms(handle) + 0.5);
    nbonds = std::intptr_t(lammps_extract_global(handle, "nbonds"));
    nangles = std::intptr_t(lammps_extract_global(handle, "nangles"));
}

int LammpsObject::write_data(int selector)
{
    std::string fname;
    if (selector == 0)
    {
        fname = "write_data " + prefixFolderOut + "/Si_results.dat";
    }
    else if (selector == 1)
    {
        fname = "write_data " + prefixFolderOut + "/Si2O3_results.dat";
    }
    else if (selector == 2)
    {
        fname = "write_data " + prefixFolderOut + "/SiO2_results.dat";
    }
    else if (selector == 3)
    {
        fname = "write_data " + prefixFolderOut + "/C_results.dat";
    }
    else if (selector == 4)
    {
        fname = "write_data " + prefixFolderOut + "/BN_results.dat";
    }
    else
    {
        exit(3);
    }

    lammps_command(handle, fname.c_str());
    return 0;
}

int LammpsObject::write_restart(int selector)
{
    std::string fname;
    if (selector == 0)
    {
        fname = "write_restart " + prefixFolderOut + "/Si_restart.restart";
    }
    else if (selector == 1)
    {
        fname = "write_restart " + prefixFolderOut + "/Si2O3_restart.restart";
    }
    else if (selector == 2)
    {
        fname = "write_restart " + prefixFolderOut + "/SiO2_restart.restart";
    }
    else if (selector == 3)
    {
        fname = "write_restart " + prefixFolderOut + "/C_restart.restart";
    }
    else if (selector == 4)
    {
        fname = "write_restart " + prefixFolderOut + "/BN_restart.restart";
    }
    else
    {
        exit(3);
    }

    lammps_command(handle, fname.c_str());

    return 0;
}

int LammpsObject::finaliseLammpsObject(int selector)
{
    std::string fname;
    if (selector == 0)
    {
        fname = "write_data " + prefixFolderOut + "/Si_results.dat";
    }
    else if (selector == 1)
    {
        fname = "write_data " + prefixFolderOut + "/Si2O3_results.dat";
    }
    else if (selector == 2)
    {
        fname = "write_data " + prefixFolderOut + "/SiO2_results.dat";
    }
    else if (selector == 3)
    {
        fname = "write_data " + prefixFolderOut + "/C_results.dat";
    }
    else if (selector == 4)
    {
        fname = "write_data " + prefixFolderOut + "/BN_results.dat";
    }

    else
    {
        exit(3);
    }

    lammps_command(handle, fname.c_str());
    lammps_close(handle);
    lammps_mpi_finalize();
    return 0;
}

void LammpsObject::runInput(std::string fname)
{
    lammps_file(handle, fname.c_str());
    nbonds = std::intptr_t(lammps_extract_global(handle, "nbonds"));
    nangles = std::intptr_t(lammps_extract_global(handle, "nangles"));
}

void LammpsObject::getatominfo(int dim)
{
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

int LammpsObject::getnAtoms()
{
    return natoms;
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
    int nbonds0, nbonds1, nangles0, nangles1;

    // unpck parameters
    int a, b, c, d, e, f;
    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];

    a += 1;
    b += 1;
    c += 1;
    d += 1;
    e += 1;
    f += 1;

    breakBond(a, e, 1, logger);
    breakBond(b, d, 1, logger);
    formBond(a, d, 1, logger);
    formBond(b, e, 1, logger);

    int e1, e11, d1, d11;
    if (networkA.nodes[e - 1].netCnxs[0] == b - 1)
    {
        e1 = networkA.nodes[e - 1].netCnxs[1] + 1;
        e11 = networkA.nodes[e - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[e - 1].netCnxs[1] == b - 1)
    {
        e1 = networkA.nodes[e - 1].netCnxs[0] + 1;
        e11 = networkA.nodes[e - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[e - 1].netCnxs[2] == b - 1)
    {
        e1 = networkA.nodes[e - 1].netCnxs[0] + 1;
        e11 = networkA.nodes[e - 1].netCnxs[1] + 1;
    }
    else
    {
        std::cout << "e not connected to b" << std::endl;
        sleep(10);
    }
    if (networkA.nodes[d - 1].netCnxs[0] == a - 1)
    {
        d1 = networkA.nodes[d - 1].netCnxs[1] + 1;
        d11 = networkA.nodes[d - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[d - 1].netCnxs[1] == a - 1)
    {
        d1 = networkA.nodes[d - 1].netCnxs[0] + 1;
        d11 = networkA.nodes[d - 1].netCnxs[2] + 1;
    }
    else if (networkA.nodes[d - 1].netCnxs[2] == a - 1)
    {
        d1 = networkA.nodes[d - 1].netCnxs[0] + 1;
        d11 = networkA.nodes[d - 1].netCnxs[1] + 1;
    }
    else
    {
        std::cout << "d not connected to a" << std::endl;
        sleep(10);
    }

    breakAngle(c, a, e);
    breakAngle(b, a, e);
    breakAngle(d, b, f);
    breakAngle(d, b, a);

    breakAngle(e1, e, a);
    breakAngle(e11, e, a);
    breakAngle(d1, d, b);
    breakAngle(d11, d, b);

    formAngle(c, a, d);
    formAngle(b, a, d);
    formAngle(f, b, e);
    formAngle(a, b, e);

    formAngle(e1, e, b);
    formAngle(e11, e, b);
    formAngle(d1, d, a);
    formAngle(d11, d, a);
}

void LammpsObject::revertGraphene(VecF<int> switchIdsA, Network networkA, LoggerPtr logger)
{

    // unpck parameters
    int a, b, c, d, e, f;
    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];

    a += 1;
    b += 1;
    c += 1;
    d += 1;
    e += 1;
    f += 1;

    breakBond(a, d, 1, logger);
    breakBond(b, e, 1, logger);
    formBond(a, e, 1, logger);
    formBond(b, d, 1, logger);

    int e1, e11, d1, d11;
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
        std::cout << "e not connected to a" << std::endl;
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
        std::cout << "d not connected to b" << std::endl;
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

void LammpsObject::switchTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB,
                                      VecF<int> switchIdsT, Network networkT, LoggerPtr logger)
{
    // unpck parameters
    int a, b, c, d, e, f;
    int u, v, w, x;
    int alpha, beta, gamma, delta, eta;
    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];
    u = switchIdsB[0];
    v = switchIdsB[1];
    w = switchIdsB[2];
    x = switchIdsB[3];

    alpha = switchIdsT[6];
    beta = switchIdsT[7];
    gamma = switchIdsT[8];
    delta = switchIdsT[9];
    eta = switchIdsT[10];

    a += 1;
    b += 1;

    alpha += 1;
    beta += 1;
    gamma += 1;
    delta += 1;
    eta += 1;

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

void LammpsObject::revertTriangleRaft(VecF<int> switchIdsA, VecF<int> switchIdsB,
                                      VecF<int> switchIdsT, LoggerPtr logger)
{
    // unpck parameters
    int a, b, c, d, e, f;
    int u, v, w, x;
    int alpha, beta, gamma, delta, eta;
    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];
    u = switchIdsB[0];
    v = switchIdsB[1];
    w = switchIdsB[2];
    x = switchIdsB[3];

    alpha = switchIdsT[6];
    beta = switchIdsT[7];
    gamma = switchIdsT[8];
    delta = switchIdsT[9];
    eta = switchIdsT[10];

    a += 1;
    b += 1;

    alpha += 1;
    beta += 1;
    gamma += 1;
    delta += 1;
    eta += 1;

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

void LammpsObject::switchBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB,
                                 VecF<int> switchIdsT, LoggerPtr logger)
{

    bool verbose = false;
    // unpck parameters
    int a, b, c, d, e, f;
    int u, v, w, x;
    int alpha, beta, gamma, delta, eta;

    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];
    u = switchIdsB[0];
    v = switchIdsB[1];
    w = switchIdsB[2];
    x = switchIdsB[3];

    alpha = switchIdsT[6];
    beta = switchIdsT[7];
    gamma = switchIdsT[8];
    delta = switchIdsT[9];
    eta = switchIdsT[10];

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

    int nSi, nSi2, nO;
    nSi = (int)(round(natoms / 3) + 0.5);
    nSi2 = (int)(round(nSi / 2) + 0.5);
    nO = natoms - nSi;

    int a_prime, a_prime_prime, b_prime, b_prime_prime;
    int a_prime_prime_prime, b_prime_prime_prime;

    int beta_prime, beta_prime_prime, gamma_prime, gamma_prime_prime;
    int eta_prime, eta_prime_prime, delta_prime, delta_prime_prime;

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

void LammpsObject::revertBilayer(VecF<int> switchIdsA, VecF<int> switchIdsB,
                                 VecF<int> switchIdsT, LoggerPtr logger)
{

    bool verbose = false;
    // unpck parameters
    int a, b, c, d, e, f;
    int u, v, w, x;
    int alpha, beta, gamma, delta, eta;

    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];
    u = switchIdsB[0];
    v = switchIdsB[1];
    w = switchIdsB[2];
    x = switchIdsB[3];

    alpha = switchIdsT[6];
    beta = switchIdsT[7];
    gamma = switchIdsT[8];
    delta = switchIdsT[9];
    eta = switchIdsT[10];

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

    int nSi, nSi2, nO;

    nSi = (int)(round(natoms / 3) + 0.5);
    nSi2 = (int)(round(nSi / 2) + 0.5);
    nO = natoms - nSi;

    int a_prime, a_prime_prime, b_prime, b_prime_prime;
    int a_prime_prime_prime, b_prime_prime_prime;

    int beta_prime, beta_prime_prime, gamma_prime, gamma_prime_prime;
    int eta_prime, eta_prime_prime, delta_prime, delta_prime_prime;

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

    std::string command;

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

void LammpsObject::switchBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT)
{
    // unpck parameters
    int a, b, c, d, e, f;
    int u, v, w, x;
    int alpha, beta, gamma, delta, eta;

    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];
    u = switchIdsB[0];
    v = switchIdsB[1];
    w = switchIdsB[2];
    x = switchIdsB[3];

    alpha = switchIdsT[6];
    beta = switchIdsT[7];
    gamma = switchIdsT[8];
    delta = switchIdsT[9];
    eta = switchIdsT[10];

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

    int nSi, nO;
    nSi = round(natoms / 3);
    nO = round(2 * natoms / 3);

    int a_prime, a_prime_prime, b_prime, b_prime_prime;
    int a_prime_prime_prime, b_prime_prime_prime;

    int beta_prime, beta_prime_prime, gamma_prime, gamma_prime_prime;
    int eta_prime, eta_prime_prime, delta_prime, delta_prime_prime;

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

void LammpsObject::revertBonds(VecF<int> switchIdsA, VecF<int> switchIdsB, VecF<int> switchIdsT)
{
    // unpck parameters
    int a, b, c, d, e, f;
    int u, v, w, x;
    int alpha, beta, gamma, delta, eta;
    a = switchIdsA[0];
    b = switchIdsA[1];
    c = switchIdsA[2];
    d = switchIdsA[3];
    e = switchIdsA[4];
    f = switchIdsA[5];
    u = switchIdsB[0];
    v = switchIdsB[1];
    w = switchIdsB[2];
    x = switchIdsB[3];

    alpha = switchIdsT[6];
    beta = switchIdsT[7];
    gamma = switchIdsT[8];
    delta = switchIdsT[9];
    eta = switchIdsT[10];
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