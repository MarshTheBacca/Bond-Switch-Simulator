
#ifndef INPUT_DATA_TPP
#define INPUT_DATA_TPP
#include "input_data.h"

/**
 * @brief Reads a value from a string and converts it to the correct type
 * @param word The word to be converted
 * @param variable The variable to be assigned
 * @param section The section of the input file
 * @tparam T The type of the value to be assigned
 */
template <typename T>
void InputData::readWord(const std::string &word, T &variable, const std::string &section) const {

    // Different conversions based on type of variable
    if constexpr (std::is_same_v<T, int>)
        try {
            variable = std::stoi(word);
        } catch (std::invalid_argument &) {
            throw std::runtime_error("Invalid integer: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    else if constexpr (std::is_same_v<T, double>) {
        try {
            variable = std::stod(word);
        } catch (std::invalid_argument &) {
            throw std::runtime_error("Invalid double: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    } else if constexpr (std::is_same_v<T, bool>) {
        try {
            variable = stringToBool(word);
        } catch (std::invalid_argument &) {
            throw std::runtime_error("Invalid boolean: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    } else if constexpr (std::is_same_v<T, StructureType>) {
        if (word == "SimpleGraphene") {
            variable = StructureType::GRAPHENE;
        } else if (word == "Silicene") {
            variable = StructureType::SILICENE;
        } else if (word == "TriangleRaft") {
            variable = StructureType::TRIANGLE_RAFT;
        } else if (word == "Bilayer") {
            variable = StructureType::BILAYER;
        } else if (word == "BoronNitride") {
            variable = StructureType::BORON_NITRIDE;
        } else {
            throw std::runtime_error("Invalid structure type: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    } else if constexpr (std::is_same_v<T, SelectionType>) {
        if (word == "Random") {
            variable = SelectionType::RANDOM;
        } else if (word == "Weighted") {
            variable = SelectionType::EXPONENTIAL_DECAY;
        } else {
            throw std::runtime_error("Invalid selection type: " + word + " in section: " + section + " on line: " + std::to_string(lineNumber));
        }
    } else if constexpr (std::is_same_v<T, std::string>) {
        variable = word;
    } else {
        throw std::invalid_argument("Cannot read word for type T. T must be int, double, bool, StructureType, SelectionType or std::string.");
    }
}

/**
 * @brief Reads a section of the input file
 * @param inputFile The input file stream
 * @param args The variables to be assigned
 * @tparam Args The types of the values to be read
 */
template <typename... Args>
void InputData::readSection(const std::string &section, Args &...args) {
    // Skip the '----' line and section title line
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    lineNumber++;
    inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    lineNumber++;

    std::string firstWord;

    // Fold expression to read data
    // For every argument in args, read a value from the input file and assign it to the argument
    ((firstWord = getFirstWord(inputFile),
      readWord(firstWord, args, section),
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
void InputData::checkInRange(const T &value, const T &lower, const T &upper, const std::string &errorMessage) const {
    if constexpr (std::is_same_v<T, int> || std::is_same_v<T, double>) {
        if (value < lower || value > upper) {
            throw std::runtime_error(errorMessage);
        }
    }
    throw std::invalid_argument("Cannot check range for type T. T must be int or double.");
}

#endif // INPUT_DATA_TPP