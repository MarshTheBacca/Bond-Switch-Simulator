// Created by olwhi on 24/07/2023, edited by Marshall Hunt 28/01/2024.
#include "lammps_object.h"
#include <filesystem>

std::string LAMMPS_FILES_PATH = std::filesystem::path("./input_files") / "lammps_files";

/**
 * @brief Default constructor for a blank Lammps Object
 */
LammpsObject::LammpsObject() = default;

/**
 * @brief Constructor for a Lammps Object
 * @param selector The selector for the structure to be created
 * @param inputFolder The folder containing the input files
 * @param loggerArg The logger object
 */
LammpsObject::LammpsObject(const LoggerPtr &loggerArg) : logger(loggerArg) {
    logger->debug("Creating Lammps Object");
    const char *lmpargv[] = {"liblammps", "-screen", "none"};
    int lmpargc = sizeof(lmpargv) / sizeof(const char *);
    logger->debug("lmpargc -> {} ", lmpargc);
    handle = lammps_open_no_mpi(lmpargc, const_cast<char **>(lmpargv), nullptr);

    if (handle == nullptr) {
        lammps_mpi_finalize();
        throw std::runtime_error("LAMMPS initialization failed");
    }
    logger->debug("Created LAMMPS handle");
    version = lammps_version(handle);
    logger->debug("LAMMPS Version: {}", version);

    std::string inputFilePath = std::filesystem::path(LAMMPS_FILES_PATH) / "lammps_script.txt";
    logger->debug("Executing LAMMPS Script: {}", inputFilePath);
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
    logger->debug("LAMMPS #nodes: {} #bonds: {} #angles: {}", natoms, nbonds, nangles);
}

/**
 * @brief Exports the network to a file
*/
void LammpsObject::writeData() {
    std::string savePath = std::filesystem::path("./output_files") / "lammps_network_result.txt";
    std::string command = "write_data " + savePath;
    lammps_command(handle, command.c_str());
}

/**
 * @brief Breaks a bond in the network
 * @param atom1 The first atom in the bond
 * @param atom2 The second atom in the bond
 * @param type The type of the bond
 * @throws std::runtime_error if the bond count doesn't decrease
*/
void LammpsObject::breakBond(const int &atom1, const int &atom2, const int &type) {
    // logger->debug("Breaking bond between {} and {} of type {}", atom1, atom2, type);
    double initialBondCount = lammps_get_thermo(handle, "bonds");

    std::string command = "group switch id " + std::to_string(atom1) + " " + std::to_string(atom2);
    lammps_command(handle, command.c_str());

    // Delete the bond
    command = "delete_bonds switch bond " + std::to_string(type) + " remove";
    lammps_command(handle, command.c_str());

    lammps_command(handle, "group switch delete");

    double finalBondCount = lammps_get_thermo(handle, "bonds");
    if (finalBondCount != initialBondCount - 1) {
        std::ostringstream oss;
        oss << "Error in Bond Counts while breaking bond between " << atom1 << " and " << atom2
            << ", initial: " << finalBondCount << " final: " << initialBondCount;
        throw std::runtime_error(oss.str());
    }
}

/**
 * @brief Forms a bond in the network
 * @param atom1 The first atom in the bond
 * @param atom2 The second atom in the bond
 * @param type The type of the bond
 * @throws std::runtime_error if the bond count doesn't increase
*/
void LammpsObject::formBond(const int &atom1, const int &atom2, const int &type) {
    // logger->debug("Forming bond between {} and {} of type {}", atom1, atom2, type);
    double initialBondCount = lammps_get_thermo(handle, "bonds");
    std::string command;
    command = "create_bonds single/bond " + std::to_string(type) + " " + std::to_string(atom1) + " " + std::to_string(atom2);
    lammps_command(handle, command.c_str());
    double finalBondCount = lammps_get_thermo(handle, "bonds");
    if (finalBondCount != initialBondCount + 1) {
        std::ostringstream oss;
        oss << "Error in Bond Counts while forming bond between " << atom1 << " and " << atom2
            << ", initial: " << finalBondCount << " final: " << initialBondCount;
        throw std::runtime_error(oss.str());
    }
}

/**
 * @brief Breaks an angle in the lattice by trying 1-2-3 first, then 3-2-1
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 * @throws std::runtime_error if the angle count doesn't decrease
 */
void LammpsObject::breakAngle(const int &atom1, const int &atom2, const int &atom3) {
    // logger->debug("Breaking angle between {}, {}, {}", atom1, atom2, atom3);
    double initialAngleCount = lammps_get_thermo(handle, "angles");
    double finalAngleCount;
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
        finalAngleCount = lammps_get_thermo(handle, "angles");
        if (finalAngleCount == initialAngleCount - 1)
            return; // The angle was successfully broken
        else if (finalAngleCount == initialAngleCount - 3) {
            // In LAMMPS, when you define a group of atoms and delete their angles, it deletes ALL the angles involving those atoms
            // This means that when you have a triangle of atoms and delete one angle, it deletes all three angles
            // So we have to add the other two angles back in
            logger->warn("Triangle angle broken, adding the other two angles back in");
            formAngle(atom1, atom3, atom2);
            formAngle(atom2, atom1, atom3);
            // Store the other two angles in a helper vector for later use
            angleHelper = {atom1, atom3, atom2, atom2, atom1, atom3};
            return;
        }

        else if (finalAngleCount == initialAngleCount - 2) {
            // This is an error prone idea, but when you delete a second angle in a triangle, it deletes both of the angles
            // added in the previous if statement, so the angle count decreases by two.
            // That's why we have to add back in the final angle, by using a helper vector that stores the other two angles
            // when the first angle was broken.
            logger->warn("Second triangle angle broken, adding the final one back in");
            if ((atom1 == angleHelper[0] && atom2 == angleHelper[1] && atom3 == angleHelper[2]) ||
                (atom1 == angleHelper[2] && atom2 == angleHelper[1] && atom3 == angleHelper[0])) {
                // If the current angle is the first angle in the angleHelper vector, then add the second angle in the angleHelper vector
                formAngle(angleHelper[3], angleHelper[4], angleHelper[5]);
                return;
            } else if ((atom1 == angleHelper[3] && atom2 == angleHelper[4] && atom3 == angleHelper[5]) ||
                       (atom1 == angleHelper[5] && atom2 == angleHelper[4] && atom3 == angleHelper[3])) {
                // If the current angle is the second angle in the angleHelper vector, then add the first angle in the angleHelper vector
                formAngle(angleHelper[0], angleHelper[1], angleHelper[2]);
                return;
            } else {
                std::ostringstream oss;
                oss << "Second triangle angle broken, but angle is not in the angleHelper vector " << atom1 << " " << atom2
                    << " " << atom3 << "not in: " << angleHelper[0] << " " << angleHelper[1] << " " << angleHelper[2] << " "
                    << angleHelper[3] << " " << angleHelper[4] << " " << angleHelper[5];
                throw std::runtime_error(oss.str());
            }
        }

        else {
            break;
        }
    }

    // If the function hasn't returned by now, the angle couldn't be broken in either direction
    std::stringstream error;
    error << "Error in Angle Counts while trying to break angle " << atom1 << " " << atom2 << " "
          << atom3 << ", initial: " << initialAngleCount << " final: " << finalAngleCount;
    throw std::runtime_error(error.str());
}

/**
 * @brief Forms an angle in the lattice
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 */
void LammpsObject::formAngle(const int &atom1, const int &atom2, const int &atom3) {
    // logger->debug("Forming angle between {}, {}, {}", atom1, atom2, atom3);
    if (atom1 == atom2 || atom1 == atom3 || atom2 == atom3) {
        std::ostringstream oss;
        oss << "Angle has one or more members that are the same ID: " << atom1 << " " << atom2 << " " << atom3;
        throw std::runtime_error(oss.str());
    }
    double initialAngleCount = lammps_get_thermo(handle, "angles");
    std::string command;
    command = "create_bonds single/angle 1 " + std::to_string(atom1) + " " + std::to_string(atom2) + " " + std::to_string(atom3);
    lammps_command(handle, command.c_str());
    double finalAngleCount = lammps_get_thermo(handle, "angles");
    if (finalAngleCount != initialAngleCount + 1) {
        std::ostringstream oss;
        oss << "Error in Angle Counts, initial while forming angle between " << atom1 << " and " << atom2 << " and " << atom3
            << ", initial: " << finalAngleCount << " final: " << initialAngleCount;
        throw std::runtime_error(oss.str());
    }
}

/**
 * @brief perform bond switch in the lattice using zero-indexed node IDs
 * @param bondBreaks Ths IDs of the bonds to be broken (1D vector of pairs)
 * @param bondMakes The IDs of the bonds to be made (1D vector of pairs)
 * @param angleBreaks The IDs of the angles to be broken (1D vector of triples)
 * @param angleMakes The IDs of the angles to be made (1D vector of triples)
 */
void LammpsObject::switchGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                                  const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes,
                                  const std::vector<double> &rotatedCoord1, const std::vector<double> &rotatedCoord2) {
    for (int i = 0; i < bondBreaks.size(); i += 2) {
        breakBond(bondBreaks[i] + 1, bondBreaks[i + 1] + 1, 1);
    }
    for (int i = 0; i < bondMakes.size(); i += 2) {
        formBond(bondMakes[i] + 1, bondMakes[i + 1] + 1, 1);
    }
    for (int i = 0; i < angleBreaks.size(); i += 3) {
        breakAngle(angleBreaks[i] + 1, angleBreaks[i + 1] + 1, angleBreaks[i + 2] + 1);
    }
    for (int i = 0; i < angleMakes.size(); i += 3) {
        formAngle(angleMakes[i] + 1, angleMakes[i + 1] + 1, angleMakes[i + 2] + 1);
    }
    int atom1ID = bondBreaks[0] + 1;
    int atom2ID = bondBreaks[2] + 1;
    setAtomCoords(atom1ID, rotatedCoord1, 2);
    setAtomCoords(atom2ID, rotatedCoord2, 2);
}
/**
 * @brief perform the opposite of a bond switch in the lattice using zero-indexed node IDs
 * @param bondBreaks Ths IDs of the bonds that have been broken (1D vector of pairs)
 * @param bondMakes The IDs of the bonds that have been made (1D vector of pairs)
 * @param angleBreaks The IDs of the angles that have been broken (1D vector of triples)
 * @param angleMakes The IDs of the angles that have been made (1D vector of triples)
 */
void LammpsObject::revertGraphene(const std::vector<int> &bondBreaks, const std::vector<int> &bondMakes,
                                  const std::vector<int> &angleBreaks, const std::vector<int> &angleMakes) {
    for (int i = 0; i < bondMakes.size(); i += 2) {
        breakBond(bondMakes[i] + 1, bondMakes[i + 1] + 1, 1);
    }
    for (int i = 0; i < bondBreaks.size(); i += 2) {
        formBond(bondBreaks[i] + 1, bondBreaks[i + 1] + 1, 1);
    }
    for (int i = 0; i < angleMakes.size(); i += 3) {
        breakAngle(angleMakes[i] + 1, angleMakes[i + 1] + 1, angleMakes[i + 2] + 1);
    }
    for (int i = 0; i < angleBreaks.size(); i += 3) {
        formAngle(angleBreaks[i] + 1, angleBreaks[i + 1] + 1, angleBreaks[i + 2] + 1);
    }
}

void LammpsObject::setCoords(std::vector<double> &newCoords, int dim) {
    if (dim != 2 && dim != 3) {
        throw std::runtime_error("Invalid dimension");
    }
    if (newCoords.size() != dim * natoms) {
        std::ostringstream oss;
        oss << "Invalid size of newCoords, expected " << dim * natoms << " got " << newCoords.size();
        throw std::runtime_error(oss.str());
    }
    lammps_scatter_atoms(handle, "x", 1, dim, newCoords.data());
}

void LammpsObject::setAtomCoords(const int &atomID, const std::vector<double> &newCoords, const int &dim) {
    if (dim != 2 && dim != 3) {
        throw std::runtime_error("Invalid dimension");
    }
    if (newCoords.size() != dim) {
        std::ostringstream oss;
        oss << "Invalid size of newCoords, expected " << dim << " got " << newCoords.size();
        throw std::runtime_error(oss.str());
    }
    auto x = (double **)lammps_extract_atom(handle, "x");
    for (int i = 0; i < dim; i++) {
        x[atomID - 1][i] = newCoords[i];
    }
}

/**
 * @brief Warns the user if the movie file already exists and starts the movie
*/
void LammpsObject::startMovie() {
    std::string movieFilePath = (std::filesystem::path("./output_files") / "simulation_movie.mpg").string();
    if (std::filesystem::exists(movieFilePath)){
        logger->warn("Movies file already exists! Overwriting simulation_movie.mpg");
    }
    std::string command = "dump myMovie all movie 1 " + movieFilePath + " type type";
    lammps_command(handle, command.c_str());
}

void LammpsObject::writeMovie() {
    lammps_command(handle, "run 0 post no");
}

void LammpsObject::stopMovie() {
    lammps_command(handle, "undump myMovie");
}

/**
 * @brief Minimise the potential energy of the network by moving atoms
 */
void LammpsObject::minimiseNetwork() {
    // Numbers are stopping tolerance for energy, stopping tolerance for force,
    // maximum number of iterations and maximum number of energy/force evaluations
    lammps_command(handle, "minimize 1.0e-6 0.0 1000000 10000000");
}

/**
 * @brief Get the potential energy of the network
 * @return The potential energy of the network
 */
double LammpsObject::getPotentialEnergy() {
    return lammps_get_thermo(handle, "pe");
}

/**
 * @brief Get the coordinates of the atoms in the network
 * @param dim the number of dimensions you want to receieve, 2 or 3
 * @return A 1D vector containing the coordinates of the atoms
 */
std::vector<double> LammpsObject::getCoords(const int &dim) const {
    if (dim != 2 && dim != 3) {
        throw std::runtime_error("Invalid dimension");
    }
    // Get the coordinates of the atoms
    std::vector<double> coords(dim * natoms);
    lammps_gather_atoms(handle, "x", 1, dim, coords.data());
    return coords;
}

/**
 * @brief Gets all the angles in the system
 * @return A 1D vector containing all the angles in the system in the form [a1atom1, a1atom2, a1atom3, a2atom1, a2atom2, a2atom3, ...]
 */
std::vector<int> LammpsObject::getAngles() const {
    const int *nangles_ptr = static_cast<int *>(lammps_extract_global(handle, "nangles"));
    int numAngles = *nangles_ptr;
    char *command_result = lammps_command(handle, "compute myAngles all property/local aatom1 aatom2 aatom3");
    auto compute_output = (double **)lammps_extract_compute(handle, "myAngles", 2, 2);

    std::vector<int> angles(numAngles * 3);
    for (int i = 0; i < numAngles; ++i) {
        for (int j = 0; j < 3; ++j) {
            angles[i * 3 + j] = static_cast<int>(compute_output[i][j]);
        }
    }
    lammps_command(handle, "uncompute myAngles");
    return angles;
}

/**
 * @brief Logs the first numLines angles in the network
 * @param numLines The number of angles to log
*/
void LammpsObject::showAngles(const int &numLines) const {
    std::vector<int> angles = getAngles();
    for (int i = 0; i < numLines; ++i) {
        logger->info("Angle: {:03} {:03} {:03}", angles[i * 3], angles[i * 3 + 1], angles[i * 3 + 2]);
    }
}
/**
 * @brief Check if an angle is not already in the system
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 * @return True if the angle is unique, false otherwise
 */
bool LammpsObject::checkAngleUnique(const int &atom1, const int &atom2, const int &atom3) const {
    std::vector<int> angles = getAngles();
    for (int i = 0; i < angles.size(); i += 3) {
        if ((angles[i] == atom1 && angles[i + 1] == atom2 && angles[i + 2] == atom3) ||
            (angles[i] == atom3 && angles[i + 1] == atom2 && angles[i + 2] == atom1)) {
            return false;
        }
    }
    return true;
}
