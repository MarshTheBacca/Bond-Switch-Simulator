#include "output_file.h"
#include <iomanip>
#include <ctime>
#include <sstream>

OutputFile::OutputFile() : file("a.out", std::ios::out | std::ios::trunc)
{
    if (!file)
    {
        throw std::runtime_error("Unable to open file: a.out");
    }
    initVariables();
}

OutputFile::OutputFile(const std::string &name) : file(name, std::ios::out | std::ios::trunc)
{
    if (!file)
    {
        std::ostringstream oss;
        oss << "Unable to open file: " << name;
        throw std::runtime_error(oss.str());
    }
    initVariables();
}

void OutputFile::initVariables(int precision, int indentSize, int sepSize, int spaceSize)
{
    file << std::fixed << std::showpoint << std::setprecision(precision);
    indent = std::string(indentSize, ' ');
    dashed = std::string(sepSize, '-');
    spacing = spaceSize;
    currIndent = 0;
}

void OutputFile::datetime(std::string message)
{
    std::time_t now = std::time(0);
    char buffer[100];
    std::strftime(buffer, 100, "%Y-%m-%d %H:%M:%S", std::localtime(&now));
    std::string time(buffer);
    file << message << time << '\n';
}

void OutputFile::separator()
{
    file << dashed << '\n';
}

LogFile::LogFile(const std::string &name) : OutputFile(name)
{
    tStart = std::chrono::high_resolution_clock::now();
    t0 = tStart;
}

LogFile::~LogFile()
{
    if (file.is_open())
    {
        tEnd = std::chrono::high_resolution_clock::now();
        double dt = std::chrono::duration_cast<std::chrono::seconds>(tEnd - tStart).count();
        dt /= 60.0;
        file << "Duration: " << dt << " minutes" << '\n';
    }
}

void LogFile::criticalError(std::string message)
{
    if (file.is_open())
    {
        file << "Critical error: " << message << '\n';
    }
    throw std::runtime_error(message);
}

double LogFile::timeElapsed()
{
    t1 = std::chrono::high_resolution_clock::now();
    double dt = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    t0 = t1;
    return dt / 1000.0;
}