#include "output_file.h"

/**
 * @brief Constructor that wipes file if it exists and creates file if it doesn't exist
 * @throw std::runtime_error if file cannot be opened
 */
OutputFile::OutputFile(const std::string &path) : file(path, std::ios::out | std::ios::trunc) {
    if (!file.is_open()) {
        std::string error_message = "Unable to open file: " + path;
        throw std::runtime_error(error_message);
    }
}

/**
 * @brief Constructor that wipes file if it exists and creates file if it doesn't exist and sets the spacing
 * @param spaceArg The spacing to be used when writing vectors
 * @throw std::runtime_error if file cannot be opened
 */
OutputFile::OutputFile(const std::string &path, const int &spaceArg) : file(path, std::ios::out | std::ios::trunc), spacing(spaceArg) {
    if (!file.is_open()) {
        std::string error_message = "Unable to open file: " + path;
        throw std::runtime_error(error_message);
    }
}

/**
 * @brief Writes the current date and time with a new line
 */
void OutputFile::writeDatetime() {
    auto now = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm localtime;
    localtime_r(&now_time_t, &localtime);
    file << std::put_time(&localtime, "%Y-%m-%d %H:%M:%S") << '\n';
}

/**
 * @brief Writes the current date and time followed by a message with a new line
 * @param message The message to be written
 */
void OutputFile::writeDatetime(const std::string &message) {
    auto now = std::chrono::system_clock::now();
    auto now_time_t = std::chrono::system_clock::to_time_t(now);
    std::tm localtime;
    localtime_r(&now_time_t, &localtime);
    file << std::put_time(&localtime, "%Y-%m-%d %H:%M:%S");
    file << " " << message << '\n';
}

/**
 * @brief Writes a string to the output file with a new line
 * @param string The string to be written
 */
void OutputFile::writeLine(const std::string &string) {
    file << string << "\n";
}

/**
 * @brief Writes a string to the output file
 * @param string The string to be written
 */
void OutputFile::write(const std::string &string) {
    file << string;
}
