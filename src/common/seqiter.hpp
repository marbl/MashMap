#pragma once

#include <string>
#include <functional>
#include "gzstream.h"

namespace seqiter {

void for_each_seq_in_file(
    const std::string& filename,
    const std::function<void(const std::string&, const std::string&)>& func) {
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
        std::cerr << "[mashmap::for_each_seq_in_file] unknown file format given to seqiter" << std::endl;
        assert(false);
        exit(1);
    }
    if (input_is_fasta) {
        while (in.good()) {
            std::string name = line.substr(1, line.find(" ")-1);
            std::string seq;
            while (std::getline(in, line)) {
                if (line[0] == '>') {
                    // this is the header of the next sequence
                    break;
                } else {
                    seq.append(line);
                }
            }
            func(name, seq);
        }
    } else if (input_is_fastq) {
        while (in.good()) {
            std::string name = line.substr(1, line.find(" ")-1);
            std::string seq;
            std::getline(in, seq); // sequence
            std::getline(in, line); // delimiter
            std::getline(in, line); // quality
            std::getline(in, line); // next header
            func(name, seq);
        }
    }
}

}
