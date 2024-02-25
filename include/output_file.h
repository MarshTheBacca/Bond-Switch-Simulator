#ifndef OUTPUT_FILE_H
#define OUTPUT_FILE_H

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

struct OutputFile {
    std::ofstream file;
    int spacing = 20;

    // Constructors
    explicit OutputFile(const std::string &name);
    explicit OutputFile(const std::string &path, const int &spaceArg);

    // Member functions
    void write(const std::string &string);
    void writeLine(const std::string &string);
    void writeDatetime();
    void writeDatetime(const std::string &message);

    // Template functions
    template <typename T>
    struct is_vector : std::false_type {};

    template <typename T>
    struct is_vector<std::vector<T>> : std::true_type {};

    template <typename T>
    void writeValue(const T &value);

    template <typename T>
    void writeVector(const std::vector<T> &vec);

    template <typename T>
    void write(const T &value);

    template <typename... Args>
    void writeValues(const Args &...args);
};

#include "output_file.tpp"
#endif // OUTPUT_FILE_H
