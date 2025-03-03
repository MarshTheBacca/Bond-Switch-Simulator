#include "file_tools.h"
#include "types.h"
#include <filesystem>
#include <fstream>
#include <unordered_set>

/**
 * Reads fixed ring IDs from a file
 * @param filePath Path to the file
 * @param logger Logger object
 * @note The file format is one integer per line
 */
std::unordered_set<int> readFixedRings(const std::filesystem::path &filePath,
                                       const LoggerPtr &logger) {
  std::ifstream fixedRingsFile(filePath);
  if (!fixedRingsFile.is_open()) {
    logger->warn(
        "Failed to open fixed rings file: {}, no fixed rings will be read",
        filePath.string());
    return {};
  }
  std::unordered_set<int> fixedRings;
  std::string line;
  while (std::getline(fixedRingsFile, line)) {
    try {
      int ringId = std::stoi(line);
      fixedRings.insert(ringId);
    } catch (const std::exception &e) {
      logger->warn("Ignoring invalid entry in file {}: {}", filePath.string(),
                   line);
    }
  }
  return fixedRings;
}