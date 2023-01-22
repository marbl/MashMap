#pragma once

#include <string>
#include <mutex>
#include <cstdio>
#include <dirent.h>

namespace yeet {

namespace temp_file {

// We use this to make the API thread-safe
std::recursive_mutex monitor;

std::string temp_dir;

bool keep_temp = false;

std::string get_dir() {
    std::lock_guard<recursive_mutex> lock(monitor);

    // Get the default temp dir from environment variables.
    if (temp_dir.empty()) {
        const char* system_temp_dir = nullptr;
        for(const char* var_name : {"TMPDIR", "TMP", "TEMP", "TEMPDIR", "USERPROFILE"}) {
            if (system_temp_dir == nullptr) {
                system_temp_dir = getenv(var_name);
            }
        }
        temp_dir = (system_temp_dir == nullptr ? "/tmp" : system_temp_dir);
    }

    return temp_dir;
}

/// Because the names are in a static object, we can delete them when
/// std::exit() is called.
struct Handler {
    std::set<std::string> filenames;
    std::string parent_directory;
    ~Handler() {
        if (!keep_temp) {
            // No need to lock in static destructor
            for (auto& filename : filenames) {
                std::remove(filename.c_str());
            }
            if (!parent_directory.empty()) {
                // There may be extraneous files in the directory still (like .fai files)
                auto directory = opendir(parent_directory.c_str());
                
                dirent* dp;
                while ((dp = readdir(directory)) != nullptr) {
                    // For every item still in it, delete it.
                    // TODO: Maybe eventually recursively delete?
                    std::remove((parent_directory + "/" + dp->d_name).c_str());
                }
                closedir(directory);
                
                // Delete the directory itself
                std::remove(parent_directory.c_str());
            }
        }
    }
} handler;

std::string create(const std::string& base,
                   const std::string& suffix) {
    std::lock_guard<recursive_mutex> lock(monitor);

    /*
    if (handler.parent_directory.empty()) {
        // Make a parent directory for our temp files
        string tmpdirname_cpp = get_dir(); // + "/wfmash-XXXXXX";
        char* tmpdirname = new char[tmpdirname_cpp.length() + 1];
        strcpy(tmpdirname, tmpdirname_cpp.c_str());
        auto got = mkdtemp(tmpdirname);
        if (got != nullptr) {
            // Save the directory we got
            handler.parent_directory = got;
        } else {
            cerr << "[wfmash]: couldn't create temp directory: " << tmpdirname << endl;
            exit(1);
        }
        delete [] tmpdirname;
    }
    */

    std::string tmpname = get_dir() + "/" + base + "XXXXXX"; // + suffix;
    // hack to use mkstemp to get us a safe temporary file name
    int fd = mkstemp(&tmpname[0]);
    if(fd != -1) {
        // we don't leave it open; we are assumed to open it again externally
        close(fd);
    } else {
        cerr << "[wfmash]: couldn't create temp file on base "
             << base << " : " << tmpname << endl;
        exit(1);
    }
    handler.filenames.insert(tmpname);
    return tmpname;
}

std::string create() {
    // No need to lock as we call this thing that locks
    return create("wfmash-", ".paf");
}

void remove(const std::string& filename) {
    std::lock_guard<recursive_mutex> lock(monitor);
    
    std::remove(filename.c_str());
    handler.filenames.erase(filename);
}

void set_dir(const std::string& new_temp_dir) {
    std::lock_guard<recursive_mutex> lock(monitor);
    
    temp_dir = new_temp_dir;
}

void set_keep_temp(bool setting) {
    keep_temp = setting;
}

}

}
