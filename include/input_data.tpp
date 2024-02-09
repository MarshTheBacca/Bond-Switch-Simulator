
#ifndef INPUT_DATA_TPP
#define INPUT_DATA_TPP
#include "input_data.h"

/**
 * @brief Reads a value from a string and converts it to the correct type
 * @param word The word to be converted
 * @param value The value to be assigned
 * @param section The section of the input file
 * @tparam T The type of the value to be assigned
 */
template <typename T>
void InputData::readValue(const std::string &word, T &value, const std::string &section, const LoggerPtr &logger) {
    // Check to see if T is of type int, double, bool, or std::string
    static_assert(std::is_same_v<T, int> || std::is_same_v<T, double> || std::is_same_v<T, bool> || std::is_same_v<T, std::string>,
                  "Invalid type for T. T must be int, double, bool, or std::string.");

    // Different conversions based on type of value
    if constexpr (std::is_same_v<T, int>)
        try {
            value = std::stoi(word);
        } catch (std::invalid_argument &) {
            logger->critical("Invalid integer: {} in section: {} on line: {}", word, section, lineNumber);
            throw std::runtime_error("Invalid integer: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    else if constexpr (std::is_same_v<T, double>) {
        try {
            value = std::stod(word);
        } catch (std::invalid_argument &) {
            logger->critical("Invalid double: {} in section: {} on line: {}", word, section, lineNumber);
            throw std::runtime_error("Invalid double: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    } else if constexpr (std::is_same_v<T, bool>) {
        try {
            value = stringToBool(word);
        } catch (std::invalid_argument &) {
            logger->critical("Invalid boolean: {} in section: {} on line: {}", word, section, lineNumber);
            throw std::runtime_error("Invalid boolean: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    } else {
        value = word;
    }
}

/**
 * @brief Reads a section of the input file
 * @param inputFile The input file stream
 * @param args The values to be read
 * @tparam Args The types of the values to be read
 */
template <typename... Args>
void InputData::readSection(std::ifstream &inputFile, const std::string &section,
                            const LoggerPtr &logger, Args &...args) {
    // Skip the '----' line and section title line
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    lineNumber++;
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    lineNumber++;

    // Initialise string stream and first word
    std::istringstream iss;
    std::string firstWord;

    // Fold expression to read data
    // For every argument in args, read a value from the input file and assign it to the argument
    ((firstWord = getFirstWord(inputFile, iss),
      readValue(firstWord, args, section, logger),
      lineNumber++),
     ...);
}

/**
 * @brief Checks if a value is in a range
 * @param value The value to be checked
 * @param lower The lower bound of the range
 * @param upper The upper bound of the range
 * @param errorMessage The error message to be thrown if the value is not in the range
 * @tparam T The type of the value to be checked
 * @throws std::runtime_error if the value is not in the range
 */
template <typename T>
void InputData::checkInRange(const T value, const T lower, const T upper, const std::string &errorMessage)
    const {
    static_assert(std::is_same_v<T, int> || std::is_same_v<T, double>,
                  "Invalid type for T. T must be int or double.");
    if (value < lower || value > upper) {
        throw std::runtime_error(errorMessage);
    }
}

#endif // INPUT_DATA_TPP