#include "input_data.h"
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
#include <type_traits>
#include <utility>

/**
 * @brief Converts a string to a boolean
 * @param str The string to be converted
 * @return The boolean value
 * @throws std::invalid_argument if the string is not "true" or "false"
 */
bool InputData::stringToBool(const std::string &str) const {
    if (str == "true")
        return true;
    else if (str == "false")
        return false;
    throw std::invalid_argument("Invalid boolean: " + str);
}

/**
 * @brief Gets the first word from a line in the input file
 * @param inputFile The input file stream
 * @param iss The input string stream
 * @return The first word from the line
 */
std::string InputData::getFirstWord(std::ifstream &inputFile,
                                    std::istringstream &iss) {
    std::string line;
    getline(inputFile, line);
    lineNumber++;
    iss.str(line);
    iss.clear();
    std::string firstWord;
    iss >> firstWord;
    return firstWord;
}

/**
 * @brief Reads the IO section of the input file
 * @param inputFile The input file stream
 * @param logger The log file
 */
void InputData::readIO(std::ifstream &inputFile, const LoggerPtr &logger) {
    readSection(inputFile, "IO", logger, outputFolder, outputFilePrefix,
                inputFolder, inputFilePrefix, isFromScratchEnabled,
                isRestartUsingLAMMPSObjectsEnabled);
}

void InputData::readNetworkProperties(std::ifstream &inputFile,
                                      const LoggerPtr &logger) {
    readSection(inputFile, "Network Properties", logger, numRings, minRingSize,
                maxRingSize, minCoordination, maxCoordination, isFixRingsEnabled,
                fixedRingsFile);
}

void InputData::readNetworkMinimisationProtocols(std::ifstream &inputFile,
                                                 const LoggerPtr &logger) {
    readSection(inputFile, "Network Minimisation Protocols", logger,
                isOpenMPIEnabled, isSimpleGrapheneEnabled, isTriangleRaftEnabled,
                isBilayerEnabled, isTersoffGrapheneEnabled, isBNEnabled,
                selectedMinimisationProtocol);
}

void InputData::readMonteCarloProcess(std::ifstream &inputFile,
                                      const LoggerPtr &logger) {
    readSection(inputFile, "Monte Carlo Process", logger, moveType, randomSeed, randomOrWeighted, weightedDecay);
}

void InputData::readMonteCarloEnergySearch(std::ifstream &inputFile,
                                           const LoggerPtr &logger) {
    readSection(inputFile, "Monte Carlo Energy Search", logger, startTemperature,
                endTemperature, temperatureIncrement, thermalisationTemperature,
                stepsPerTemperature, initialThermalisationSteps);
}

void InputData::readPotentialModel(std::ifstream &inputFile,
                                   const LoggerPtr &logger) {
    readSection(inputFile, "Potential Model", logger, maximumBondLength, maximumAngle);
}

void InputData::readAnalysis(std::ifstream &inputFile,
                             const LoggerPtr &logger) {
    readSection(inputFile, "Analysis", logger, analysisWriteFrequency,
                isWriteSamplingStructuresEnabled, structureWriteFrequency);
}

void InputData::readOutput(std::ifstream &inputFile, const LoggerPtr &logger) {
    readSection(inputFile, "Output", logger, ljPairsCalculationDistance);
}

void InputData::checkInSet(const std::string &value,
                           const std::set<std::string, std::less<>> &validValues,
                           const std::string &errorMessage) const {
    if (validValues.count(value) == 0) {
        throw std::runtime_error(errorMessage);
    }
}

void InputData::checkFileExists(const std::string &filename) const {
    std::ifstream file(filename);
    if (!file) {
        throw std::runtime_error("File does not exist: " + filename);
    }
}

void InputData::validate() {
    // Network Properties
    checkInRange(numRings, 1, INT_MAX, "Number of rings must be at least 1");
    checkInRange(minRingSize, 3, INT_MAX, "Minimum ring size must be at least 3");
    checkInRange(maxRingSize, minRingSize, INT_MAX, "Maximum ring size must be at least the minimum ring size");
    checkInRange(minCoordination, 1, INT_MAX, "Minimum coordination must be at least 1");
    checkInRange(maxCoordination, minCoordination, INT_MAX, "Maximum coordination must be at least the minimum coordination");
    if (isFixRingsEnabled) {
        checkFileExists(fixedRingsFile);
    }

    // Minimisation Protocols
    std::map<int, std::pair<bool *, std::string>> protocolMap = {
        {1, {&isSimpleGrapheneEnabled, "Simple Graphene"}},
        {2, {&isTriangleRaftEnabled, "Triangle Raft"}},
        {3, {&isTriangleRaftEnabled, "Bilayer"}},
        {4, {&isTersoffGrapheneEnabled, "Tersoff Graphene"}},
        {5, {&isBNEnabled, "BN"}}};

    if (protocolMap.count(selectedMinimisationProtocol) == 1) {
        auto &[isEnabled, protocolName] = protocolMap[selectedMinimisationProtocol];
        if (!*isEnabled) {
            throw std::runtime_error("Selected minimisation protocol is " + std::to_string(selectedMinimisationProtocol) +
                                     " but " + protocolName + " is disabled");
        }
    } else {
        throw std::runtime_error("Selected minimisation protocol, " + std::to_string(selectedMinimisationProtocol) +
                                 " is out of range");
    }

    // Monte Carlo Process
    checkInSet(moveType, {"switch", "mix"}, "Invalid move type: " + moveType + " must be either 'switch' or 'mix'");
    checkInRange(randomSeed, 0, INT_MAX, "Random seed must be at least 0");
    checkInSet(randomOrWeighted, {"random", "weighted"}, "Invalid random or weighted: " + randomOrWeighted + " must be either 'random' or 'weighted'");

    // Analysis
    checkInRange(analysisWriteFrequency, 0, 1000, "Analysis write frequency must be between 0 and 1000");

    // Output
    checkInRange(ljPairsCalculationDistance, 0, INT_MAX, "LJ pairs calculation distance must be at least 0");
}

/**
 * @brief Reads the input file
 * @param filePath The path to the input file
 * @param logger The log file
 */
InputData::InputData(const std::string &filePath, const LoggerPtr &logger) {
    // Open the input file
    std::ifstream inputFile(filePath);

    // Check if the file was opened successfully
    if (!inputFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + filePath);
    }
    logger->info("Reading input file: " + filePath);

    // Skip the title line
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    lineNumber++;

    // Read sections
    readIO(inputFile, logger);
    readNetworkProperties(inputFile, logger);
    readNetworkMinimisationProtocols(inputFile, logger);
    readMonteCarloProcess(inputFile, logger);
    readMonteCarloEnergySearch(inputFile, logger);
    readPotentialModel(inputFile, logger);
    readAnalysis(inputFile, logger);

    // Validate input data
    logger->info("Validating input data...");
    validate();
}