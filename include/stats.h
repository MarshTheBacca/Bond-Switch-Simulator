#ifndef D33FA5A5_564D_4EFF_9BD8_9D4AD9108DC5
#define D33FA5A5_564D_4EFF_9BD8_9D4AD9108DC5

class Stats {
public:
  Stats() = default;

  void incrementSwitches() { switches++; }
  void incrementAcceptedSwitches() { acceptedSwitches++; }
  void incrementFailedBondLengthChecks() { failedBondLengthChecks++; }
  void incrementFailedAngleChecks() { failedAngleChecks++; }
  void incrementFailedEnergyChecks() { failedEnergyChecks++; }

  int getSwitches() const { return switches; }
  int getAcceptedSwitches() const { return acceptedSwitches; }
  int getFailedBondLengthChecks() const { return failedBondLengthChecks; }
  int getFailedAngleChecks() const { return failedAngleChecks; }
  int getFailedEnergyChecks() const { return failedEnergyChecks; }
  void reset() {
    switches = 0;
    acceptedSwitches = 0;
    failedBondLengthChecks = 0;
    failedAngleChecks = 0;
    failedEnergyChecks = 0;
  }

private:
  int switches = 0;               // Number of switches performed
  int acceptedSwitches = 0;       // Number of switches accepted
  int failedBondLengthChecks = 0; // Number of failed bond length checks
  int failedAngleChecks = 0;      // Number of failed angle checks
  int failedEnergyChecks = 0;     // Number of failed energy checks
};

#endif /* D33FA5A5_564D_4EFF_9BD8_9D4AD9108DC5 */
