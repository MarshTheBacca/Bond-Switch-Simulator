#ifndef OUTPUT_FILE_H
#define OUTPUT_FILE_H

#include <fstream>
#include <string>
#include <chrono>
#include <vector>
#include "vecf.h"

class OutputFile
{
public:
    int currIndent;

    explicit OutputFile(const std::string &name);
    virtual ~OutputFile() = default;

    void initVariables(int precision = 8, int indentSize = 4, int sepSize = 60, int spaceSize = 20);
    void datetime(const std::string &message = "");
    void separator();
    template <typename T>
    void write(T val);
    template <typename T, typename U>
    void write(T val0, U val1);
    template <typename T>
    void writeRowVector(T vec);

protected:
    std::ofstream file;
    std::string indent;
    std::string dashed;
    int spacing;
};

#include "output_file.tpp"
#endif // OUTPUT_FILE_H
