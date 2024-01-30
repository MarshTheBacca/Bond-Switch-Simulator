#ifndef INPUT_DATA_H
#define INPUT_DATA_H

#include <string>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <limits>
#include <type_traits>
#include <set>
#include "output_file.h"
#include <spdlog/spdlog.h>

using LoggerPtr = std::shared_ptr<spdlog::logger>;

struct InputData
{
    // Used for error messages
    int lineNumber = 0;

    // IO Data
    std::string outputFolder;
    std::string outputFilePrefix;
    std::string inputFolder;
    std::string inputFilePrefix;
    bool isFromScratchEnabled;
    bool isRestartUsingLAMMPSObjectsEnabled;

    // Network Properties Data
    int numRings;
    int minRingSize;
    int maxRingSize;
    int minCoordination;
    int maxCoordination;
    bool isFixRingsEnabled;
    std::string fixedRingsFile;

    // Minimisation Protocols Data
    bool isOpenMPIEnabled;
    bool isSimpleGrapheneEnabled;
    bool isTriangleRaftEnabled;
    bool isBilayerEnabled;
    bool isTersoffGrapheneEnabled;
    bool isBNEnabled;
    int selectedMinimisationProtocol;

    // Monte Carlo Process Data
    std::string moveType;
    int randomSeed;
    bool isSpiralEnabled;
    int spiralRadius;
    std::string randomOrWeighted;

    // Monte Carlo Energy Search Data
    double startTemperature;
    double endTemperature;
    double temperatureIncrement;
    double thermalisationTemperature;
    int stepsPerTemperature;
    int initialThermalisationSteps;

    // Potential Model Data
    double harmonicBondForceConstant;
    double harmonicAngleForceConstant;
    double harmonicGeometryConstraint;
    bool isMaintainConvexityEnabled;

    // Geometry Optimisation Data
    int monteCarloLocalMaxIterations;
    int globalMinimisationMaxIterations;
    double tauBacktrackingParameter;
    double tolerance;
    int localRegionSize;

    // Analysis Data
    int analysisWriteFrequency;
    bool isWriteSamplingStructuresEnabled;
    int structureWriteFrequency;

    // Output Data
    int ljPairsCalculationDistance;

    // Declare template functions
    template <typename T>
    void readValue(const std::string &word, T &value, const std::string &section, const LoggerPtr &logger);
    template <typename... Args>
    void readSection(std::ifstream &inputFile, const std::string &section, const LoggerPtr &logger, Args &...args);
    template <typename T>
    void checkInRange(const T value, const T lower, const T upper, const std::string &errorMessage);

    // Declare non-template functions
    bool stringToBool(const std::string &str);
    std::string getFirstWord(std::ifstream &inputFile, std::istringstream &iss);
    void readIO(std::ifstream &inputFile, const LoggerPtr &logger);
    void readNetworkProperties(std::ifstream &inputFile, const LoggerPtr &logger);
    void readNetworkMinimisationProtocols(std::ifstream &inputFile, const LoggerPtr &logger);
    void readMonteCarloProcess(std::ifstream &inputFile, const LoggerPtr &logger);
    void readMonteCarloEnergySearch(std::ifstream &inputFile, const LoggerPtr &logger);
    void readPotentialModel(std::ifstream &inputFile, const LoggerPtr &logger);
    void readGeometryOptimisation(std::ifstream &inputFile, const LoggerPtr &logger);
    void readAnalysis(std::ifstream &inputFile, const LoggerPtr &logger);
    void readOutput(std::ifstream &inputFile, const LoggerPtr &logger);
    void checkInSet(const std::string &value, const std::set<std::string> &validValues, const std::string &errorMessage);
    void checkFileExists(const std::string &filename);
    void validate();
    InputData(const std::string &filePath, const LoggerPtr &logger);
};

#include "input_data.tpp"
#endif // INPUT_DATA_H
