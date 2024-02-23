#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include "output_file.h"
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

using LoggerPtr = std::shared_ptr<spdlog::logger>;

enum class StructureType { GRAPHENE,
                           SILICENE,
                           TRIANGLE_RAFT,
                           BILAYER,
                           BORON_NITRIDE };

enum class SelectionType {
    RANDOM,
    EXPONENTIAL_DECAY
};

struct InputData {
    // Used for error messages
    int lineNumber = 0;
    std::ifstream inputFile;
    // IO Data
    std::string outputFolder;
    std::string outputFilePrefix;
    std::string inputFolder;
    std::string inputFilePrefix;

    bool isFromScratchEnabled;

    // Network Properties Data
    int numRings;
    int minRingSize;
    int maxRingSize;
    int minCoordination;
    int maxCoordination;
    bool isFixRingsEnabled;

    // Minimisation Protocols Data
    bool isOpenMPIEnabled;
    StructureType structureType;

    // Monte Carlo Process Data
    int randomSeed;
    SelectionType randomOrWeighted;
    double weightedDecay;

    // Monte Carlo Energy Search Data
    double startTemperature;
    double endTemperature;
    double temperatureIncrement;
    double thermalisationTemperature;
    int stepsPerTemperature;
    int initialThermalisationSteps;

    // Potential Model Data
    double maximumBondLength;
    double maximumAngle;

    // Analysis Data
    int analysisWriteInterval;
    bool writeMovie;

    LoggerPtr logger;

    // Declare template functions
    template <typename T>
    void readWord(const std::string &word, T &variable, const std::string &section) const;
    template <typename... Args>
    void readSection(const std::string &section, Args &...args);
    template <typename T>
    void checkInRange(const T &value, const T &lower, const T &upper, const std::string &errorMessage) const;

    // Declare non-template functions
    bool stringToBool(const std::string &str) const;
    std::string getFirstWord();

    void readIO();
    void readNetworkProperties();
    void readNetworkMinimisationProtocols();
    void readMonteCarloProcess();
    void readMonteCarloEnergySearch();
    void readPotentialModel();
    void readAnalysis();

    void checkFileExists(const std::string &filename) const;
    void validate() const;
    InputData(const std::string &filePath, const LoggerPtr &logger);
};

#include "input_data.tpp"
#endif // INPUT_DATA_H
