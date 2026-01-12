#ifndef OUTPUT_FILE_H
#define OUTPUT_FILE_H

#include "stats.h"
#include <chrono>
#include <ctime>
#include <fstream>
#include <map>
#include <string>
#include <type_traits>
#include <vector>

template <typename T> struct is_vector : std::false_type {};

template <typename T> struct is_vector<std::vector<T>> : std::true_type {};

template <typename T> struct is_map : std::false_type {};

template <typename K, typename V>
struct is_map<std::map<K, V>> : std::true_type {};

template <typename T> struct is_nested_map : std::false_type {};

template <typename K, typename VK, typename VV>
struct is_nested_map<std::map<K, std::map<VK, VV>>> : std::true_type {};

struct OutputFile {
  std::ofstream file;
  int spacing = 20;

  // Constructors
  explicit OutputFile(const std::string &name);
  explicit OutputFile(const std::string &path, const int spaceArg);
  // Move constructor
  OutputFile(OutputFile &&other) noexcept
      : file(std::move(other.file)), spacing(other.spacing) {
    // Leave other in a valid state
    other.spacing = 0;
  }
  // Move assignment operator
  OutputFile &operator=(OutputFile &&other) noexcept {
    if (this != &other) {
      if (file.is_open()) {
        file.close();
      }
      file = std::move(other.file);
      spacing = other.spacing;
      other.spacing = 0;
    }
    return *this;
  }
  // Destructor
  ~OutputFile() {
    if (file.is_open()) {
      file.close();
    }
  }

  // Member functions
  void write(const std::string &string);
  void writeLine(const std::string &string);
  void writeDatetime();
  void writeDatetime(const std::string &message);
  void writeFooter(
      const Stats &stats, const bool networkConsistent,
      const std::chrono::time_point<std::chrono::high_resolution_clock> &start);

  // Template functions
  template <typename T> void writeValue(const T &value) { file << value; }

  template <typename T> void writeVector(const std::vector<T> &vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
      file << vec[i];
      if (i < vec.size() - 1) {
        file << ';';
      }
    }
  }

  template <typename K, typename V> void writeMap(const std::map<K, V> &map) {
    for (auto it = map.begin(); it != map.end(); ++it) {
      file << it->first << ':' << it->second;
      if (std::next(it) != map.end()) {
        file << ";";
      }
    }
  }

  template <typename K, typename VK, typename VV>
  void writeNestedMap(const std::map<K, std::map<VK, VV>> &map) {
    for (auto it = map.begin(); it != map.end(); ++it) {
      file << it->first << '!';
      writeMap(it->second);
      if (std::next(it) != map.end()) {
        file << "&";
      }
    }
  }

  template <typename T> void write(const T &value) {
    if constexpr (is_vector<T>::value) {
      writeVector(value);
    } else if constexpr (is_nested_map<T>::value) {
      writeNestedMap(value);
    } else if constexpr (is_map<T>::value) {
      writeMap(value);
    } else {
      writeValue(value);
    }
  }

  template <typename... Args> void writeValues(const Args &...args) {
    int n = 0;
    ((write(args), file << (n++ < sizeof...(Args) - 1 ? ',' : '\n')), ...);
  }
};

#endif // OUTPUT_FILE_H
