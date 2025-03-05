#ifndef SWITCH_MOVE_H
#define SWITCH_MOVE_H

#include "types.h"
#include <array>
#include <cstdint>
#include <exception>
#include <string>
#include <unordered_set>

// custom exception
class SwitchMoveException : public std::exception {
public:
  explicit SwitchMoveException(const std::string &message) : message(message) {}
  const char *what() const noexcept override { return message.c_str(); }

private:
  std::string message;
};

struct SwitchMove {
  const std::array<uint16_t, 2> selectedBaseBond;
  const std::array<uint16_t, 2> selectedRingBond;
  const std::array<std::array<uint16_t, 2>, 2> bondBreaks;
  const std::array<std::array<uint16_t, 2>, 2> bondMakes;
  const std::array<std::array<uint16_t, 3>, 8> angleBreaks;
  const std::array<std::array<uint16_t, 3>, 8> angleMakes;
  const std::array<std::array<uint16_t, 2>, 2> ringBondBreakMake;

  const std::unordered_set<uint16_t> involvedBaseNodes;

  // default constructor
  SwitchMove();

  SwitchMove(const std::array<uint16_t, 2> selectedBaseBond,
             const std::array<uint16_t, 2> selectedRingBond,
             const std::array<std::array<uint16_t, 2>, 2> &bondBreaks,
             const std::array<std::array<uint16_t, 2>, 2> &bondMakes,
             const std::array<std::array<uint16_t, 3>, 8> &angleBreaks,
             const std::array<std::array<uint16_t, 3>, 8> &angleMakes,
             const std::array<std::array<uint16_t, 2>, 2> &ringBondBreakMake,
             const std::unordered_set<uint16_t> &involvedBaseNodes)
      : selectedBaseBond(selectedBaseBond), selectedRingBond(selectedRingBond),
        bondBreaks(bondBreaks), bondMakes(bondMakes), angleBreaks(angleBreaks),
        angleMakes(angleMakes), ringBondBreakMake(ringBondBreakMake),
        involvedBaseNodes(involvedBaseNodes) {}
};

#endif // SWITCH_MOVE_H