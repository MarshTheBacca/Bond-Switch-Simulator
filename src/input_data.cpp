#include "input_data.h"

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
    getline(inputFile, line);
    lineNumber++;
    std::istringstream iss(line);
    std::string firstWord;
    iss >> firstWord;
    return firstWord;
}

/**
 * @brief Reads the IO section of the input file
 * @param inputFile The input file stream
 * @param logger The log file
 */
void InputData::readIO() {
    readSection("IO", outputFolder, outputFilePrefix,
                inputFolder, inputFilePrefix, isFromScratchEnabled);
}

void InputData::readNetworkProperties() {
    readSection("Network Properties", numRings, minRingSize,
                maxRingSize, minCoordination, maxCoordination, isFixRingsEnabled);
}

void InputData::readNetworkMinimisationProtocols() {
    readSection("Network Minimisation Protocols",
                isOpenMPIEnabled, structureType);
}

void InputData::readMonteCarloProcess() {
    readSection("Monte Carlo Process", randomSeed, randomOrWeighted, weightedDecay);
}

void InputData::readMonteCarloEnergySearch() {
    readSection("Monte Carlo Energy Search", startTemperature,
                endTemperature, temperatureIncrement, thermalisationTemperature,
                stepsPerTemperature, initialThermalisationSteps);
}

void InputData::readPotentialModel() {
    readSection("Potential Model", maximumBondLength, maximumAngle);
}

void InputData::readAnalysis() {
    readSection("Analysis", analysisWriteFrequency, writeMovie);
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
    checkInRange(numRings, 1, INT_MAX, "Number of rings must be at least 1");
    checkInRange(minRingSize, 3, INT_MAX, "Minimum ring size must be at least 3");
    checkInRange(maxRingSize, minRingSize, INT_MAX, "Maximum ring size must be at least the minimum ring size");
    checkInRange(minCoordination, 1, INT_MAX, "Minimum coordination must be at least 1");
    checkInRange(maxCoordination, minCoordination, INT_MAX, "Maximum coordination must be at least the minimum coordination");
    if (isFixRingsEnabled) {
        checkFileExists(inputFolder + "/fixed_rings.dat");
    }
    // Monte Carlo Process
    checkInRange(randomSeed, 0, INT_MAX, "Random seed must be at least 0");

    // Energy search
    if (endTemperature < startTemperature && temperatureIncrement > 0) {
        throw std::runtime_error("Temperature increment must be negative if end temperature is less than start temperature");
    } else if (endTemperature > startTemperature && temperatureIncrement < 0) {
        throw std::runtime_error("Temperature increment must be positive if end temperature is greater than start temperature");
    }

    // Analysis
    checkInRange(analysisWriteFrequency, 0, 1000, "Analysis write frequency must be between 0 and 1000");
}

/**
 * @brief Reads the input file
 * @param filePath The path to the input file
 * @param logger The log file
 */
InputData::InputData(const std::string &filePath, const LoggerPtr &logger) : inputFile(filePath) {

    // Check if the file was opened successfully
    if (!inputFile.is_open()) {
        throw std::runtime_error("Unable to open file: " + filePath);
    }
    logger->info("Reading input file: " + filePath);

    // Skip the title line
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    lineNumber++;

    // Read sections
    readIO();
    readNetworkProperties();
    readNetworkMinimisationProtocols();
    readMonteCarloProcess();
    readMonteCarloEnergySearch();
    readPotentialModel();
    readAnalysis();

    // Validate input data
    logger->info("Validating input data...");
    validate();
}