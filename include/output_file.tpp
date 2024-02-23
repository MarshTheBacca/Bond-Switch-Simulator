#ifndef OUTPUT_FILE_TPP
#define OUTPUT_FILE_TPP

#include "output_file.h"


/**
 * @brief Writes two values to the output file with a new line
 *
 *
 */
template <typename... Args>
void OutputFile::writeValues(Args... args) {
    (file << ... << args) << '\n';
}

/**
 * @brief Writes a vector to the output file with a new line
 * @param vec The vector to be written
 * @tparam T The type of the vector to be written
 */
template <typename T>
void OutputFile::writeVector(const std::vector<T> &vec) {
    for (int i = 0; i < vec.size(); ++i) {
        file << std::setw(spacing) << std::left << vec[i];
    }
    file << '\n';
}

#endif // OUTPUT_FILE_TPP