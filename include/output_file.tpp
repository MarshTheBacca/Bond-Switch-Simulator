#ifndef OUTPUT_FILE_TPP
#define OUTPUT_FILE_TPP

#include <iomanip>
#include <ctime>
#include <sstream>
#include "output_file.h"

/**
 * @brief Writes a value to the output file with a new line
 * @param val The value to be written
 * @tparam T The type of the value to be written
 */
template <typename T>
void OutputFile::write(T val)
{
    file << std::string(currIndent, ' ') << val << '\n';
}

/**
 * @brief Writes two values to the output file with a new line
 * @param val0 The first value to be written
 * @param val1 The second value to be written
 * @tparam T The type of the first value to be written
 * @tparam U The type of the second value to be written
 */
template <typename T, typename U>
void OutputFile::write(T val0, U val1)
{
    file << std::string(currIndent, ' ') << val0 << " " << val1 << '\n';
}

/**
 * @brief Writes a row vector to the output file with a new line
 * @param vec The vector to be written
 * @tparam T The type of the vector to be written
 */
template <typename T>
void OutputFile::writeRowVector(T vec)
{
    // static_assert(std::is_same<decltype(vec.n), int>::value, "Type T must have an integer member 'n'");
    // static_assert(std::is_same<decltype(vec[0]), double>::value || std::is_same<decltype(vec[0]), int>::value,
    //            "Vector must have a subscript operator returning double or int");
    for (int i = 0; i < vec.n; ++i)
    {
        file << std::setw(spacing) << std::left << vec[i];
    }
    file << '\n';
}

#endif // OUTPUT_FILE_TPP