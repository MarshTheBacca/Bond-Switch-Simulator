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

    OutputFile();
    explicit OutputFile(const std::string &name);
    virtual ~OutputFile() = default;

    void initVariables(int precision = 8, int indentSize = 4, int sepSize = 60, int spaceSize = 20);
    void datetime(std::string message = "");
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

class LogFile : public OutputFile
{
public:
    explicit LogFile(const std::string &name);
    ~LogFile();

    void criticalError(std::string message);
    double timeElapsed();

private:
    std::chrono::high_resolution_clock::time_point tStart;
    std::chrono::high_resolution_clock::time_point t0;
    std::chrono::high_resolution_clock::time_point t1;
    std::chrono::high_resolution_clock::time_point tEnd;
};

#include "output_file.tpp"
#endif // OUTPUT_FILE_H