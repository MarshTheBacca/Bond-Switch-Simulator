#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include "output_file.h"
#include "types.h"
#include <climits>
#include <fstream>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <spdlog/spdlog.h>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

enum class StructureType {
  GRAPHENE,
  SILICENE,
  TRIANGLE_RAFT,
  BILAYER,
  BORON_NITRIDE
};

enum class SelectionType { RANDOM, EXPONENTIAL_DECAY };

struct InputData {
  // Used for error messages
  int lineNumber = 0;
  std::ifstream inputFile;

  // Network Restrictions Data
  int minRingSize;
  int maxRingSize;
  double maximumBondLength;
  double maximumAngle;
  bool isFixRingsEnabled;

  // Bond Selection Process Data
  int randomSeed;
  SelectionType randomOrWeighted;
  double weightedDecay;

  // Temperature Schedule Data
  double thermalisationTemperature;
  double annealingStartTemperature;
  double annealingEndTemperature;
  int thermalisationSteps;
  int annealingSteps;

  // Analysis Data
  int analysisWriteInterval;
  bool writeMovie;

  LoggerPtr logger;

  InputData(const std::string &filePath, const LoggerPtr &logger);

  // Declare template functions
  template <typename T>
  void readWord(const std::string &word, T &variable,
                const std::string &section) const;
  template <typename... Args>
  void readSection(const std::string &section, Args &...args);
  template <typename T>
  void checkInRange(const T &value, const T &lower, const T &upper,
                    const std::string &errorMessage) const;

  // Declare non-template functions
  bool stringToBool(const std::string &str) const;
  std::string getFirstWord();

  void readNetworkRestrictions();
  void readBondSelectionProcess();
  void readTemperatureSchedule();
  void readAnalysis();

  void checkFileExists(const std::string &filename) const;
  void validate() const;
};

#include "input_data.tpp"
#endif // INPUT_DATA_H
