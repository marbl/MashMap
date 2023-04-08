#pragma once

#include <sys/stat.h>

namespace fs {
/**
 * Check if a file exists
 * @return true if and only if the file exists, false else
 */
bool file_exists(const std::string& file) {
    struct stat buf;
    return (stat(file.c_str(), &buf) == 0);
}

}
