#include "output_file.h"
#include <chrono>

/**
 * @brief Constructor that wipes file if it exists and creates file if it
 * doesn't exist
 * @throw std::runtime_error if file cannot be opened
 */
OutputFile::OutputFile(const std::string &path)
    : file(path, std::ios::out | std::ios::trunc) {
  if (!file.is_open()) {
    std::string error_message = "Unable to open file: " + path;
    throw std::runtime_error(error_message);
  }
}

/**
 * @brief Constructor that wipes file if it exists and creates file if it
 * doesn't exist and sets the spacing
 * @param spaceArg The spacing to be used when writing vectors
 * @throw std::runtime_error if file cannot be opened
 */
OutputFile::OutputFile(const std::string &path, const int spaceArg)
    : file(path, std::ios::out | std::ios::trunc), spacing(spaceArg) {
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
void OutputFile::write(const std::string &string) { file << string; }

void OutputFile::writeFooter(
    const Stats &stats, const bool networkConsistent,
    const std::chrono::time_point<std::chrono::high_resolution_clock> &start) {
  auto end = std::chrono::high_resolution_clock::now();
  auto duration =
      std::chrono::duration_cast<std::chrono::milliseconds>(end - start) /
      1000.0;
  this->writeLine(
      "The following line is a few statistics about the simulation");
  this->writeLine(
      "Number of attempted switches, Number of accepted switches, Number of "
      "failed angle checks, Number of failed bond length checks, Number of "
      "failed energy checks, Monte Carlo acceptance, Total run time (s), "
      "Average time per step (us), Network Consistent");
  this->writeValues(
      stats.getSwitches(), stats.getAcceptedSwitches(),
      stats.getFailedAngleChecks(), stats.getFailedBondLengthChecks(),
      stats.getFailedEnergyChecks(),
      (double)stats.getAcceptedSwitches() / stats.getSwitches(),
      duration.count(), duration.count() / stats.getSwitches() * 1000.0,
      networkConsistent ? "true" : "false");
}