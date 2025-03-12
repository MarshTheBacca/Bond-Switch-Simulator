#ifndef LAMMPS_OBJECT_H
#define LAMMPS_OBJECT_H
#include "network.h"
#include "switch_move.h"
#include "types.h"
#include <vector>

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

  std::vector<double> getCoords(const int dim) const;
  void setCoords(std::vector<double> &newCoords, int dim);
  void setAtomCoords(const int atomID, const std::vector<double> &newCoords,
                     const int dim);

  void breakBond(const int atom1, const int atom2, const int type);
  void formBond(const int atom1, const int atom2, const int type);
  void breakAngle(const int atom1, const int atom2, const int atom3);
  void formAngle(const int atom1, const int atom2, const int atom3);
  void breakAndFormAngles(int atom3, int atom1, int atom5, int atom4, int e1,
                          int e11, int d1, int d11);

  void switchGraphene(const SwitchMove &switchMove,
                      const std::vector<double> &rotatedCoord1,
                      const std::vector<double> &rotatedCoord2);

  void revertGraphene(const SwitchMove &switchMove);
  std::vector<int> getAngles() const;

  void writeData();
  void startMovie();
  void writeMovie();
  void stopMovie();

  void showAngles(const int numLines) const;
  bool checkAngleUnique(const int atom1, const int atom2,
                        const int atom3) const;
};

#endif // LAMMPS_OBJECT_H
