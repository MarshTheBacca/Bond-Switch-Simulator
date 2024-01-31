#include "output_file.h"
#include <iomanip>
#include <ctime>
#include <sstream>

OutputFile::OutputFile(const std::string &name) : file(name, std::ios::out | std::ios::trunc)
{
    if (!file.is_open())
    {
        std::string error_message = "Unable to open file: " + name;
        throw std::runtime_error(error_message);
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

void OutputFile::datetime(const std::string &message)
{
    std::time_t now = std::time(0);
    struct tm timeinfo;
    localtime_r(&now, &timeinfo);
    char buffer[100];
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", &timeinfo);
    std::string time(buffer);
    file << message << time << '\n';
}

void OutputFile::separator()
{
    file << dashed << '\n';
}