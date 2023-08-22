#pragma once

#include <string>
#include <functional>
#include <cassert>
#include <unordered_set>
#include "gzstream.h"
#ifdef USE_HTSLIB
#include <htslib/faidx.h>
#endif

namespace seqiter {

bool fai_index_exists(const std::string& filename) {
  // Check if .fai file exists
  std::ifstream f(filename + ".fai");
  return f.good();
}

void for_each_seq_in_file(
    const std::string& filename,
    const std::unordered_set<std::string>& keep_seq,
    const std::string& keep_prefix,
    const std::function<void(const std::string&, const std::string&)>& func) {

#ifdef USE_HTSLIB
    if ((!keep_seq.empty() || !keep_prefix.empty())
          && fai_index_exists(filename)) {
        // Use index
        // Handle keep_prefix
        std::unordered_set<std::string> found_seq;
        auto faid = fai_load(filename.c_str());
        if (!keep_prefix.empty()) {
            // iterate over lines in the fai file
            std::string line;
            std::ifstream in(filename + ".fai");
            while (std::getline(in, line)) {
                std::string name = line.substr(0, line.find("\t"));
                char* seq;
                int64_t len = 0;
                if (strncmp(name.c_str(), keep_prefix.c_str(), keep_prefix.size()) == 0
                    || keep_seq.find(name) != keep_seq.end()) {
                    found_seq.insert(name);
                    seq = faidx_fetch_seq64(
                        faid, name.c_str(), 0, INT_MAX, &len);
                    func(name, std::string(seq, len));
                    free(seq);
                } else {
                    func(name, "");
                }
            }
        }
        // Handle keep_seq
        for (const auto& name : keep_seq) {
            if (found_seq.find(name) == found_seq.end())
            {
                std::cerr << "[wfmash::for_each_seq_in_file] could not fetch " << name << " from index" << std::endl;
            }
        }
        fai_destroy(faid); // Free FAI index
    } else 
#endif
    {
        // no index available
        // detect file type
        bool input_is_fasta = false;
        bool input_is_fastq = false;
        std::string line;
        igzstream in(filename.c_str());
        std::getline(in, line);
        if (line[0] == '>') {
            input_is_fasta = true;
        } else if (line[0] == '@') {
            input_is_fastq = true;
        } else {
            std::cerr << "[wfmash::for_each_seq_in_file] unknown file format given to seqiter" << std::endl;
            assert(false);
            exit(1);
        }
        if (input_is_fasta) {
            while (in.good()) {
                std::string name = line.substr(1, line.find(" ")-1);
                std::string seq;
                bool keep = (keep_prefix.empty() || name.substr(0, keep_prefix.length()) == keep_prefix)
                    && (keep_seq.empty() || keep_seq.find(name) != keep_seq.end());
                while (std::getline(in, line)) {
                    if (line[0] == '>') {
                        // this is the header of the next sequence
                        break;
                    } else {
                        if (keep) {
                            seq.append(line);
                        }
                    }
                }
                func(name, seq);
            }
        } else if (input_is_fastq) {
            while (in.good()) {
                std::string name = line.substr(1, line.find(" ")-1);
                std::string seq;
                bool keep = (keep_prefix.empty() || name.substr(0, keep_prefix.length()) == keep_prefix)
                    && (keep_seq.empty() || keep_seq.find(name) != keep_seq.end());
                std::getline(in, seq); // sequence
                std::getline(in, line); // delimiter
                std::getline(in, line); // quality
                std::getline(in, line); // next header
                func(name, keep ? seq : "");
            }
        }
    }
}

} // namespace seqiter
