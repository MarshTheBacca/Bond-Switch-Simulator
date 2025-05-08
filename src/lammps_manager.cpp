#include "lammps_manager.h"
#include "library.h" // LAMMPS Library
#include "switch_move.h"
#include <filesystem>
#include <format>

const std::string LAMMPS_FILES_PATH =
    std::filesystem::path("./input_files") / "lammps_files";

/**
 * @brief Default constructor for a blank LAMMPS Manager
 */
LAMMPSManager::LAMMPSManager() = default;

/**
 * @brief Constructor for a LAMMPS Manager
 * @param selector The selector for the structure to be created
 * @param inputFolder The folder containing the input files
 * @param loggerArg The logger object
 */
LAMMPSManager::LAMMPSManager(const LoggerPtr &loggerArg) : logger(loggerArg) {
  logger->debug("Creating LAMMPS Manager");
  std::array<char *, 3> lammpsArguments = {"liblammps", "-screen", "none"};
  auto numLammpsArguments = static_cast<int>(lammpsArguments.size());
  handle =
      lammps_open_no_mpi(numLammpsArguments, lammpsArguments.data(), nullptr);

  if (handle == nullptr) {
    lammps_mpi_finalize();
    throw std::runtime_error("LAMMPS initialization failed");
  }
  logger->debug("Created LAMMPS handle");
  version = lammps_version(handle);
  logger->debug("LAMMPS Version: {}", version);

  std::string inputFilePath =
      std::filesystem::path(LAMMPS_FILES_PATH) / "lammps_script.txt";
  logger->debug("Executing LAMMPS Script: {}", inputFilePath);
  lammps_file(handle, inputFilePath.c_str());
  numAtoms = (int)(lammps_get_natoms(handle) + 0.5);
  if (const auto nbonds_ptr =
          static_cast<const int *>(lammps_extract_global(handle, "nbonds"));
      nbonds_ptr) {
    nbonds = *nbonds_ptr;
  } else {
    logger->error("Failed to extract number of bonds");
  }

  if (const auto nangles_ptr =
          static_cast<int *>(lammps_extract_global(handle, "nangles"));
      nangles_ptr) {
    nangles = *nangles_ptr;
  } else {
    logger->error("Failed to extract number of angles");
  }
  logger->debug("LAMMPS #nodes: {} #bonds: {} #angles: {}", numAtoms, nbonds,
                nangles);
}

/**
 * @brief Exports the network to a file
 */
void LAMMPSManager::writeData() {
  std::string savePath =
      std::filesystem::path("./output_files") / "lammps_network_result.txt";
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
void LAMMPSManager::breakBond(const uint16_t atom1, const uint16_t atom2,
                              const int type) {
  size_t initialBondCount = this->getBondCount();
  // Define a group of atoms called `switch`
  lammps_command(handle,
                 std::format("group switch id {} {}", atom1, atom2).c_str());

  // Delete the bond described by the group `switch`
  lammps_command(
      handle, std::format("delete_bonds switch bond {} remove", type).c_str());

  // Remove the group definition
  lammps_command(handle, "group switch delete");
  size_t finalBondCount = this->getBondCount();
  if (finalBondCount != initialBondCount - 1) {
    throw std::runtime_error(std::format("Error in Bond Counts while breaking "
                                         "bond between {} and {}, initial: {} "
                                         "final: {}",
                                         atom1, atom2, finalBondCount,
                                         initialBondCount));
  }
}

/**
 * @brief Forms a bond in the network
 * @param atom1 The first atom in the bond
 * @param atom2 The second atom in the bond
 * @param type The type of the bond
 * @throws std::runtime_error if the bond count doesn't increase
 */
void LAMMPSManager::formBond(const uint16_t atom1, const uint16_t atom2,
                             const int type) {
  size_t initialBondCount = this->getBondCount();
  std::string command =
      std::format("create_bonds single/bond {} {} {}", type, atom1, atom2);
  lammps_command(handle, command.c_str());
  size_t finalBondCount = this->getBondCount();
  if (finalBondCount != initialBondCount + 1) {
    throw std::runtime_error(
        std::format("Error in Bond Counts while forming bond between {} and "
                    "{}, initial count: {} final count: {}",
                    atom1, atom2, finalBondCount, initialBondCount));
  }
}

/**
 * @brief Breaks an angle in the lattice by trying 1-2-3 first, then 3-2-1
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 * @throws std::runtime_error if the angle count doesn't decrease
 */
void LAMMPSManager::breakAngle(const uint16_t atom1, const uint16_t atom2,
                               const uint16_t atom3) {
  // Get the initial number of angles
  const size_t initialAngleCount = this->getAngleCount();
  // Define a group `switch` describing the atoms in the angle
  lammps_command(
      handle,
      std::format("group switch id {} {} {}", atom1, atom2, atom3).c_str());
  // Delete the angle described by the group `switch`
  lammps_command(handle, "delete_bonds switch angle 1 remove");

  // Delete the user defined group
  lammps_command(handle, "group switch delete");

  // Check if the number of angles has decreased
  const size_t finalAngleCount = this->getAngleCount();
  if (finalAngleCount == initialAngleCount - 1) {
    // The angle was successfully broken
    return;
  } else if (finalAngleCount == initialAngleCount - 3) {
    // In LAMMPS, when you define a group of atoms and delete their angles,
    // it deletes ALL the angles involving those atoms This means that when
    // you have a triangle of atoms and delete one angle, it deletes all
    // three angles So we have to add the other two angles back in
    logger->warn("Triangle angle broken, adding the other two angles back in");
    formAngle(atom1, atom3, atom2);
    formAngle(atom2, atom1, atom3);
    // Store the other two angles in a helper vector for later use
    angleHelper = {atom1, atom3, atom2, atom2, atom1, atom3};
    return;
  } else if (finalAngleCount == initialAngleCount - 2) {
    // This is an error prone idea, but when you delete a second angle in a
    // triangle, it deletes both of the angles added in the previous if
    // statement, so the angle count decreases by two. That's why we have to
    // add back in the final angle, by using a helper vector that stores the
    // other two angles when the first angle was broken.
    logger->warn("Second triangle angle broken, adding the final one back in");
    if ((atom1 == angleHelper[0] && atom2 == angleHelper[1] &&
         atom3 == angleHelper[2]) ||
        (atom1 == angleHelper[2] && atom2 == angleHelper[1] &&
         atom3 == angleHelper[0])) {
      // If the current angle is the first angle in the angleHelper vector,
      // then add the second angle in the angleHelper vector
      formAngle(angleHelper[3], angleHelper[4], angleHelper[5]);
      return;
    } else if ((atom1 == angleHelper[3] && atom2 == angleHelper[4] &&
                atom3 == angleHelper[5]) ||
               (atom1 == angleHelper[5] && atom2 == angleHelper[4] &&
                atom3 == angleHelper[3])) {
      // If the current angle is the second angle in the angleHelper vector,
      // then add the first angle in the angleHelper vector
      formAngle(angleHelper[0], angleHelper[1], angleHelper[2]);
      return;
    } else {
      throw std::runtime_error(std::format(
          "Second triangle angle broken, but angle {}-{}-{} is not in the "
          "angleHelper vector {}-{}-{}, {}-{}-{}",
          atom1, atom2, atom3, angleHelper[0], angleHelper[1], angleHelper[2],
          angleHelper[3], angleHelper[4], angleHelper[5]));
    }
  }

  // If the function hasn't returned by now, the angle couldn't be broken in
  // either direction
  throw std::runtime_error(std::format("Error in Angle Counts while trying to "
                                       "break angle {}-{}-{}, initial: {} "
                                       "final: {}",
                                       atom1, atom2, atom3, initialAngleCount,
                                       finalAngleCount));
}

/**
 * @brief Forms an angle in the lattice
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 */
void LAMMPSManager::formAngle(const uint16_t atom1, const uint16_t atom2,
                              const uint16_t atom3) {
  if (atom1 == atom2 || atom1 == atom3 || atom2 == atom3) {
    throw std::runtime_error(std::format("Angle has one or more members that "
                                         "are the same ID: {} {} {}",
                                         atom1, atom2, atom3));
  }
  size_t initialAngleCount = this->getAngleCount();
  lammps_command(handle, std::format("create_bonds single/angle 1 {} {} {}",
                                     atom1, atom2, atom3)
                             .c_str());
  size_t finalAngleCount = this->getAngleCount();
  if (finalAngleCount != initialAngleCount + 1) {
    throw std::runtime_error(std::format("Error in Angle Counts while forming "
                                         "angle between {} and {} and {}, "
                                         "initial: {} final: {}",
                                         atom1, atom2, atom3, finalAngleCount,
                                         initialAngleCount));
  }
}

/**
 * @brief perform bond switch in the lattice using zero-indexed node IDs
 */
void LAMMPSManager::performSwitch(const SwitchMove &switchMove,
                                  const std::array<double, 2> &rotatedCoord1,
                                  const std::array<double, 2> &rotatedCoord2) {
  std::ranges::for_each(switchMove.bondBreaks,
                        [this](const std::array<uint16_t, 2> &bond) {
                          breakBond(bond[0] + 1, bond[1] + 1, 1);
                        });
  std::ranges::for_each(switchMove.bondMakes,
                        [this](const std::array<uint16_t, 2> &bond) {
                          formBond(bond[0] + 1, bond[1] + 1, 1);
                        });
  std::ranges::for_each(switchMove.angleBreaks,
                        [this](const std::array<uint16_t, 3> &angle) {
                          breakAngle(angle[0] + 1, angle[1] + 1, angle[2] + 1);
                        });
  std::ranges::for_each(switchMove.angleMakes,
                        [this](const std::array<uint16_t, 3> &angle) {
                          formAngle(angle[0] + 1, angle[1] + 1, angle[2] + 1);
                        });

  int atom1ID = switchMove.bondBreaks[0][0] + 1;
  int atom2ID = switchMove.bondBreaks[1][0] + 1;
  setAtomCoords(atom1ID, rotatedCoord1);
  setAtomCoords(atom2ID, rotatedCoord2);
}

/**
 * @brief perform the opposite of a bond switch in the lattice using
 * zero-indexed node IDs
 */
void LAMMPSManager::revertSwitch(const SwitchMove &switchMove) {
  std::ranges::for_each(switchMove.bondMakes,
                        [this](const std::array<uint16_t, 2> &bond) {
                          breakBond(bond[0] + 1, bond[1] + 1, 1);
                        });
  std::ranges::for_each(switchMove.bondBreaks,
                        [this](const std::array<uint16_t, 2> &bond) {
                          formBond(bond[0] + 1, bond[1] + 1, 1);
                        });
  std::ranges::for_each(switchMove.angleMakes,
                        [this](const std::array<uint16_t, 3> &angle) {
                          breakAngle(angle[0] + 1, angle[1] + 1, angle[2] + 1);
                        });
  std::ranges::for_each(switchMove.angleBreaks,
                        [this](const std::array<uint16_t, 3> &angle) {
                          formAngle(angle[0] + 1, angle[1] + 1, angle[2] + 1);
                        });
}

void LAMMPSManager::setCoords(
    const std::vector<std::array<double, 2>> &newCoords) {
  if (newCoords.size() != numAtoms) {
    throw std::runtime_error(std::format("Invalid size of newCoords, expected "
                                         "{} got {}",
                                         numAtoms, newCoords.size()));
  }
  // Flatten the 2D vector into a 1D vector, which LAMMPS expects
  std::vector<double> coords1D(numAtoms * 2);
  for (size_t i = 0; i < newCoords.size(); ++i) {
    coords1D[2 * i] = newCoords[i][0];
    coords1D[2 * i + 1] = newCoords[i][1];
  }
  lammps_scatter_atoms(
      // The LAMMPS handle
      handle,
      // x means positions
      "x",
      // Atom of type 1
      1,
      // 2 values per atom (2D)
      2,
      // The coordinate values as a 1D array
      coords1D.data());
}

void LAMMPSManager::setAtomCoords(const int atomID,
                                  const std::array<double, 2> &newCoords) {
  auto x = (double **)lammps_extract_atom(handle, "x");
  for (int i = 0; i < 2; i++) {
    x[atomID - 1][i] = newCoords[i];
  }
}

/**
 * @brief Warns the user if the movie file already exists and starts the movie
 */
void LAMMPSManager::startMovie() {
  std::string movieFilePath =
      (std::filesystem::path("./output_files") / "simulation_movie.mpg")
          .string();
  if (std::filesystem::exists(movieFilePath)) {
    logger->warn(
        "Movies file already exists! Overwriting simulation_movie.mpg");
  }
  std::string command =
      "dump myMovie all movie 1 " + movieFilePath + " type type";
  lammps_command(handle, command.c_str());
}

void LAMMPSManager::writeMovie() { lammps_command(handle, "run 0 post no"); }

void LAMMPSManager::stopMovie() { lammps_command(handle, "undump myMovie"); }

/**
 * @brief Minimise the potential energy of the network by moving atoms
 */
void LAMMPSManager::minimiseNetwork() {
  // Numbers are stopping tolerance for energy, stopping tolerance for force,
  // maximum number of iterations and maximum number of energy/force
  // evaluations
  lammps_command(handle, "minimize ${etol} ${ftol} ${maxiter} ${maxeval}");
}

/**
 * @brief Get the potential energy of the network
 * @return The potential energy of the network
 */
double LAMMPSManager::getPotentialEnergy() {
  return lammps_get_thermo(handle, "pe");
}

/**
 * @brief Get the coordinates of the atoms in the network
 * @return A 1D vector containing the coordinates of the atoms
 */
std::vector<std::array<double, 2>> LAMMPSManager::getCoords() const {
  // Get the coordinates of the atoms
  std::vector<double> coords(2 * numAtoms);
  lammps_gather_atoms(
      // LAMMPS handle
      handle,
      // Get coordinates
      "x",
      // 1 means doubles, 0 means integers
      1,
      // The number of values per atom
      2,
      // The buffer to populate
      coords.data());
  if (coords.size() != numAtoms * 2) {
    throw std::runtime_error(
        std::format("LAMMPS returned coords of size {}, expected {}",
                    coords.size(), numAtoms * 2));
  }
  std::vector<std::array<double, 2>> coords2D(numAtoms);
  for (int i = 0; i < numAtoms; i++) {
    coords2D[i] = {coords[i * 2], coords[i * 2 + 1]};
  }
  return coords2D;
}

/**
 * @brief Gets all the angles in the system
 * @return A 1D vector containing all the angles in the system in the form
 * [a1atom1, a1atom2, a1atom3, a2atom1, a2atom2, a2atom3, ...]
 */
std::vector<int> LAMMPSManager::getAngles() const {
  const int *nangles_ptr =
      static_cast<int *>(lammps_extract_global(handle, "nangles"));
  int numAngles = *nangles_ptr;
  lammps_command(handle,
                 "compute myAngles all property/local aatom1 aatom2 aatom3");
  auto compute_output =
      (double **)lammps_extract_compute(handle, "myAngles", 2, 2);

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
void LAMMPSManager::showAngles(const size_t numLines) const {
  std::vector<int> angles = getAngles();
  for (size_t i = 0; i < numLines; ++i) {
    logger->info("Angle: {:03} {:03} {:03}", angles[i * 3], angles[i * 3 + 1],
                 angles[i * 3 + 2]);
  }
}
/**
 * @brief Check if an angle is not already in the system
 * @param atom1 The first atom in the angle
 * @param atom2 The second atom in the angle
 * @param atom3 The third atom in the angle
 * @return True if the angle is unique, false otherwise
 */
bool LAMMPSManager::checkAngleUnique(const uint16_t atom1, const uint16_t atom2,
                                     const uint16_t atom3) const {
  std::vector<int> angles = getAngles();
  for (int i = 0; i < angles.size(); i += 3) {
    if ((angles[i] == atom1 && angles[i + 1] == atom2 &&
         angles[i + 2] == atom3) ||
        (angles[i] == atom3 && angles[i + 1] == atom2 &&
         angles[i + 2] == atom1)) {
      return false;
    }
  }
  return true;
}

size_t LAMMPSManager::getBondCount() const {
  return static_cast<size_t>(lammps_get_thermo(handle, "bonds"));
}

size_t LAMMPSManager::getAngleCount() const {
  return static_cast<size_t>(lammps_get_thermo(handle, "angles"));
}