#ifndef OUTPUT_FILE_TPP
#define OUTPUT_FILE_TPP
#include "output_file.h"

template <typename T>
void OutputFile::writeValue(const T &value) {
    file << value;
}

template <typename T>
void OutputFile::writeVector(const std::vector<T> &vec) {
    for (size_t i = 0; i < vec.size(); ++i) {
        file << vec[i];
        if (i < vec.size() - 1) {
            file << ';';
        }
    }
}

template <typename K, typename V>
void OutputFile::writeMap(const std::map<K, V> &map) {
    for (auto it = map.begin(); it != map.end(); ++it) {
        file << it->first << ':' << it->second;
        if (std::next(it) != map.end()) {
            file << ";";
        }
    }
}

template <typename K, typename VK, typename VV>
void OutputFile::writeNestedMap(const std::map<K, std::map<VK, VV>> &map) {
    for (auto it = map.begin(); it != map.end(); ++it) {
        file << it->first << '!';
        writeMap(it->second);
        if (std::next(it) != map.end()) {
            file << "&";
        }
    }
}

template <typename T>
void OutputFile::write(const T &value) {
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

template <typename... Args>
void OutputFile::writeValues(const Args &...args) {
    int n = 0;
    ((write(args), file << (n++ < sizeof...(Args) - 1 ? ',' : '\n')), ...);
}

#endif // OUTPUT_FILE_TPP