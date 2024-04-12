#include "input_data.h"
#include <filesystem>

/**
 * @brief Reads the input file
 * @param filePath The path to the input file
 * @param loggerArg The log file
 */
InputData::InputData(const std::string &filePath, const LoggerPtr &loggerArg) : inputFile(filePath), logger(loggerArg) {

    // Check if the file was opened successfully
    if (!inputFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + filePath);
    }
    logger->debug("Reading input file: " + filePath);

    // Skip the title line
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    lineNumber++;

    // Read sections
    readNetworkRestrictions();
    readBondSelectionProcess();
    readTemperatureSchedule();
    readAnalysis();

    // Validate input data
    logger->debug("Validating input data...");
    validate();
}

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
 * @brief Gets the first word from a line in the input file and increment line number
 * @return The first word from the line as a string
 */
std::string InputData::getFirstWord() {
    std::string line;
    std::getline(inputFile, line);
    lineNumber++;
    std::istringstream iss(line);
    std::string firstWord;
    iss >> firstWord;
    return firstWord;
}

void InputData::readNetworkRestrictions() {
    readSection("Network Restrictions",  minRingSize, maxRingSize, maximumBondLength, maximumAngle, isFixRingsEnabled);
}

void InputData::readBondSelectionProcess() {
    readSection("Bond Selection Process", randomSeed, randomOrWeighted, weightedDecay);
}

void InputData::readTemperatureSchedule() {
    readSection("Temperature", thermalisationTemperature, annealingStartTemperature,
                annealingEndTemperature, thermalisationSteps, annealingSteps);
}

void InputData::readAnalysis() {
    readSection("Analysis", analysisWriteInterval, writeMovie);
}

/**
 * @brief Checks if a file exists
 * @param path The path of the file
 * @throws std::runtime_error if the file does not exist
 */
void InputData::checkFileExists(const std::string &path) const {
    std::ifstream file(path);
    if (!file) {
        throw std::runtime_error("File does not exist: " + path);
    }
}

/**
 * @brief Validates the input data
 * @throws std::runtime_error if the input data is invalid
 */
void InputData::validate() const {
    // Network Properties
    checkInRange(minRingSize, 3, INT_MAX, "Minimum ring size must be at least 3");
    checkInRange(maxRingSize, minRingSize, INT_MAX, "Maximum ring size must be at least the minimum ring size");

    if (isFixRingsEnabled) {
        checkFileExists(std::filesystem::path("./input_files") / "bss_network" / "fixed_rings.txt");
    }
    // Monte Carlo Process
    checkInRange(randomSeed, 0, INT_MAX, "Random seed must be at least 0");

    // Monte Carlo Energy Searcg
    checkInRange(annealingSteps, 0, INT_MAX, "Annealing steps must be at least 0");
    checkInRange(thermalisationSteps, 0, INT_MAX, "Thermalisation steps must be at least 0");
    checkInRange(maximumBondLength, 0.0, std::numeric_limits<double>::max(), "Maximum bond length must be at least 0");
    checkInRange(maximumAngle, 0.0, 360.0, "Maximum angle must be between 0 and 360");

    // Analysis
    if (analysisWriteInterval < 0) {
        throw std::runtime_error("Analysis write interval must be at least 0");
    }
    if (writeMovie && thermalisationSteps + annealingSteps > 2000) {
        throw std::runtime_error("Cannot write a movie file for more than 2000 steps because the file would be enormous");
    }
}
