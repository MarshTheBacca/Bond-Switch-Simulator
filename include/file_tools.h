#ifndef FILE_TOOLS_H
#define FILE_TOOLS_H

#include "types.h"
#include <cstdint>
#include <filesystem>
#include <unordered_set>

std::unordered_set<uint16_t>
readFixedRings(const std::filesystem::path &filePath, const LoggerPtr &logger);

#endif // FILE_TOOLS_H