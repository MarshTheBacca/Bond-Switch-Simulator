#ifndef TYPES_H
#define TYPES_H

#include "spdlog/spdlog.h"
#include <memory>

using LoggerPtr = std::shared_ptr<spdlog::logger>;
enum class Direction { CLOCKWISE, ANTICLOCKWISE };

#endif // TYPES_H