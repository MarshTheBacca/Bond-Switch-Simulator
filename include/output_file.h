#ifndef OUTPUT_FILE_H
#define OUTPUT_FILE_H

#include <chrono>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

struct OutputFile {
    std::ofstream file;
    int spacing = 20;

    explicit OutputFile(const std::string &name);
    explicit OutputFile(const std::string &path, const int &spaceArg);

    void write(const std::string &string);
    void writeLine(const std::string &string);

    template <typename... Args>
    void writeValues(Args... args);

    template <typename T>
    void writeVector(const std::vector<T> &vec);

    void writeDatetime();
    void writeDatetime(const std::string &message = "");
};

#include "output_file.tpp"
#endif // OUTPUT_FILE_H
