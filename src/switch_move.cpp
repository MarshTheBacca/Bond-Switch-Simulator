#include "switch_move.h"
#include "types.h"
#include <cstdint>

// Default constructor
SwitchMove::SwitchMove()
    : selectedBaseBond{}, selectedRingBond{}, bondBreaks{}, bondMakes{},
      angleBreaks{}, angleMakes{}, ringBondBreakMake{}, involvedBaseNodes{} {}
