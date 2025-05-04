#ifndef LAMMPS_MANAGER_H
#define LAMMPS_MANAGER_H
#include "network.h"
#include "switch_move.h"
#include "types.h"
#include <array>
#include <vector>

struct LAMMPSManager {
  Network networkA;
  Network networkB;

  void *handle;
  int version;
  int numAtoms = 0;
  int nbonds = 0;
  int nangles = 0;
  double *bonds = nullptr;

  std::vector<uint16_t> angleHelper = std::vector<uint16_t>(6);

  LoggerPtr logger;

  LAMMPSManager();
  explicit LAMMPSManager(const LoggerPtr &loggerArg);

  void minimiseNetwork();
  double getPotentialEnergy();

  std::vector<std::array<double, 2>> getCoords() const;
  void setCoords(const std::vector<std::array<double, 2>> &newCoords);
  void setAtomCoords(const int atomID, const std::array<double, 2> &newCoords);

  void breakBond(const uint16_t atom1, const uint16_t atom2, const int type);
  void formBond(const uint16_t atom1, const uint16_t atom2, const int type);
  void breakAngle(const uint16_t atom1, const uint16_t atom2,
                  const uint16_t atom3);
  void formAngle(const uint16_t atom1, const uint16_t atom2,
                 const uint16_t atom3);

  void performSwitch(const SwitchMove &switchMove,
                     const std::array<double, 2> &rotatedCoord1,
                     const std::array<double, 2> &rotatedCoord2);

  void revertSwitch(const SwitchMove &switchMove);
  std::vector<int> getAngles() const;

  void writeData();
  void startMovie();
  void writeMovie();
  void stopMovie();

  void showAngles(const size_t numLines) const;
  bool checkAngleUnique(const uint16_t atom1, const uint16_t atom2,
                        const uint16_t atom3) const;
  size_t getBondCount() const;
  size_t getAngleCount() const;
};

#endif // LAMMPS_MANAGER_H
