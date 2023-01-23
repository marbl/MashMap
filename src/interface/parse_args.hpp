#pragma once

#include <unistd.h>

#include "common/args.hxx"

#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"

#include "interface/temp_file.hpp"
#include "common/utils.hpp"

#include "mashmap_git_version.hpp"

// If the MASHMAP_GIT_VERSION doesn't exist at all, define a placeholder
#ifndef MASHMAP_GIT_VERSION
#define MASHMAP_GIT_VERSION "not-from-git"
#endif

namespace yeet {

struct Parameters {
    bool approx_mapping = false;
    bool remapping = false;
    //bool align_input_paf = false;
};

int64_t handy_parameter(const std::string& value) {
    auto is_a_float = [](const std::string s) {
        return !s.empty() && s.find_first_not_of("0123456789.") == std::string::npos && std::count(s.begin(), s.end(), '.') < 2;
    };

    uint64_t str_len = value.length();
    uint8_t exp = 0;
    if (value[str_len-1] == 'k' || value[str_len-1] == 'K') {
        exp = 3;
        --str_len;
    } else if (value[str_len-1] == 'm' || value[str_len-1] == 'M') {
        exp = 6;
        --str_len;
    } else if (value[str_len-1] == 'g' || value[str_len-1] == 'G') {
        exp = 9;
        --str_len;
    }

    const std::string tmp = value.substr(0, str_len);
    return is_a_float(tmp) ? (int)(stof(tmp) * pow(10, exp)) : -1;
}

void parse_args(int argc,
                char** argv,
                skch::Parameters& map_parameters,
                yeet::Parameters& yeet_parameters) {

    args::ArgumentParser parser("mashmap: a pangenome-scale mapper, " + std::string(MASHMAP_GIT_VERSION));
    parser.helpParams.width = 100;
    parser.helpParams.showTerminator = false;

    args::Group mandatory_opts(parser, "[ MANDATORY OPTIONS ]");
    args::Positional<std::string> target_sequence_file(mandatory_opts, "target", "alignment target/reference sequence file");

    args::Group io_opts(parser, "[ Files IO Options ]");
    args::PositionalList<std::string> query_sequence_files(io_opts, "queries", "query sequences file");
    args::ValueFlag<std::string> query_sequence_file_list(io_opts, "queries", "alignment queries files list", {'Q', "query-file-list"});

    args::Group mapping_opts(parser, "[ Mapping Options ]");
    args::ValueFlag<std::string> index_load_file(io_opts, "FILENAME", "Filename of index file to load", {'A', "index-load-file"});
    args::ValueFlag<std::string> index_save_file(io_opts, "FILENAME", "Filename of index file to save", {'z', "index-save-file"});
    args::ValueFlag<std::string> segment_length(mapping_opts, "N", "segment seed length for mapping [default: 5k]", {'s', "segment-length"});
    args::ValueFlag<std::string> block_length(mapping_opts, "N", "keep merged mappings supported by homologies of this total length [default: 5*segment-length]", {'l', "block-length"});
    args::ValueFlag<std::string> chain_gap(mapping_opts, "N", "chain mappings closer than this distance in query and target, retaining mappings in best chain [default: 100k]", {'c', "chain-gap"});
    args::ValueFlag<int> kmer_size(mapping_opts, "N", "kmer size [default: 19]", {'k', "kmer"});
    args::ValueFlag<float> kmer_pct_threshold(mapping_opts, "%", "ignore the top % most-frequent kmers [default: 0.001]", {'H', "kmer-threshold"});
    args::ValueFlag<uint32_t> num_mappings_for_segments(mapping_opts, "N", "number of mappings to retain for each segment [default: 1]", {'n', "num-mappings-for-segment"});
    args::ValueFlag<uint32_t> num_mappings_for_short_seq(mapping_opts, "N", "number of mappings to retain for each sequence shorter than segment length [default: 1]", {'S', "num-mappings-for-short-seq"});
    args::Flag skip_self(mapping_opts, "", "skip self mappings when the query and target name is the same (for all-vs-all mode)", {'X', "skip-self"});
    args::ValueFlag<char> skip_prefix(mapping_opts, "C", "skip mappings when the query and target have the same prefix before the last occurrence of the given character C", {'Y', "skip-prefix"});
    args::Flag no_split(mapping_opts, "no-split", "disable splitting of input sequences during mapping [default: enabled]", {'N',"no-split"});
    args::ValueFlag<float> map_pct_identity(mapping_opts, "%", "percent identity in the mashmap step [default: 90]", {'p', "map-pct-id"});
    args::Flag drop_low_map_pct_identity(mapping_opts, "K", "drop mappings with estimated identity below --map-pct-id=%", {'K', "drop-low-map-id"});
    args::Flag keep_low_align_pct_identity(mapping_opts, "A", "keep alignments with gap-compressed identity below --map-pct-id=% x 0.75", {'O', "keep-low-align-id"});
    args::Flag stage2_full_scan(mapping_opts, "stage2-full-scan", "scan full candidate regions for best minhash instead of just using the point with the highest intersection [default: disabled]", {'F',"s2-full-scan"});
    args::Flag disable_topANI_filter(mapping_opts, "disable-top-ANI-filter", "Do not use the threshold filtering for stage 1 of mapping", {'D', "disable-topANI-filter"});
    args::ValueFlag<float> map_ani_threshold(mapping_opts, "%", "ANI difference threshold for stage 1 filtering [default: 0.0]", {'T', "s1-ani-thresh"});
    args::ValueFlag<float> map_ani_threshold_conf(mapping_opts, "%", "Confidence for ANI difference threshold for stage 1 filtering [default: 0.999]", {'C', "s1-ani-thresh-conf"});
    args::Flag no_filter(mapping_opts, "MODE", "disable mapping filtering", {'f', "no-filter"});
    args::ValueFlag<double> map_sparsification(mapping_opts, "FACTOR", "keep this fraction of mappings", {'x', "sparsify-mappings"});
    args::Flag no_merge(mapping_opts, "no-merge", "don't merge consecutive segment-level mappings (NOT FULLY IMPLEMENTED)", {'M', "no-merge"});
    args::ValueFlag<int64_t> sketchSize(mapping_opts, "N", "Number of sketch elements [default 25]", {'J', "sketch-size"});
    //args::ValueFlag<std::string> path_high_frequency_kmers(mapping_opts, "FILE", " input file containing list of high frequency kmers", {'H', "high-freq-kmers"});

    args::Group general_opts(parser, "[ General Options ]");
    args::ValueFlag<std::string> tmp_base(general_opts, "PATH", "base name for temporary files [default: `pwd`]", {'B', "tmp-base"});
    args::Flag keep_temp_files(general_opts, "", "keep intermediate files", {'Z', "keep-temp"});
    //args::Flag show_progress(general_opts, "show-progress", "write alignment progress to stderr", {'P', "show-progress"});

    args::Group threading_opts(parser, "[ Threading ]");
    args::ValueFlag<int> thread_count(threading_opts, "N", "use this many threads during parallel steps", {'t', "threads"});

    args::Group program_info_opts(parser, "[ Program Information ]");
    args::Flag version(program_info_opts, "version", "show version number and github commit hash", {'v', "version"});
    args::HelpFlag help(program_info_opts, "help", "display this help menu", {'h', "help"});

    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        exit(0);
        //return; // 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        exit(1);
        //return; // 1;
    }
    if (argc==1) {
        std::cout << parser;
        exit(1);
        //return; // 1;
    }

    if (version) {
        std::cerr << MASHMAP_GIT_VERSION << std::endl;
        exit(0);
    }

    if (skip_self) {
        map_parameters.skip_self = true;
    } else {
        map_parameters.skip_self = false;
    }

    if (skip_prefix) {
        map_parameters.skip_prefix = true;
        map_parameters.prefix_delim = args::get(skip_prefix);
    } else {
        map_parameters.skip_prefix = false;
        map_parameters.prefix_delim = '\0';
    }

    if (target_sequence_file) {
        map_parameters.refSequences.push_back(args::get(target_sequence_file));
    }
    map_parameters.referenceSize = skch::CommonFunc::getReferenceSize(map_parameters.refSequences);

    if (query_sequence_files) {
        for (auto& q : args::get(query_sequence_files)) {
            map_parameters.querySequences.push_back(q);
        }
    }
    if (query_sequence_file_list) {
        skch::parseFileList(args::get(query_sequence_file_list), map_parameters.querySequences);
    }

    if (index_save_file) {
        map_parameters.saveIndexFilename = args::get(index_save_file);
    } else {
        map_parameters.saveIndexFilename = "";
    }

    if (index_load_file) {
        map_parameters.loadIndexFilename = args::get(index_load_file);
    } else {
        map_parameters.loadIndexFilename = "";
    }

    // If there are no queries, go in all-vs-all mode with the sequences specified in `target_sequence_file`
    if (target_sequence_file && map_parameters.querySequences.empty()) {
        map_parameters.skip_self = true;
        map_parameters.querySequences.push_back(map_parameters.refSequences.back());
    }

    map_parameters.alphabetSize = 4;

    if (no_filter) {
        map_parameters.filterMode = skch::filter::NONE;
    } else {
        if (map_parameters.skip_self || map_parameters.skip_prefix) {
            // before we set skch::filter::ONETOONE here
            // but this does not provide a clear benefit in all-to-all
            // as it sometimes introduces cases of over-filtering
            map_parameters.filterMode = skch::filter::MAP;
        } else {
            map_parameters.filterMode = skch::filter::MAP;
        }
    }

    if (map_sparsification) {
        if (args::get(map_sparsification) == 1) {
            // overflows
            map_parameters.sparsity_hash_threshold = std::numeric_limits<uint64_t>::max();
        } else {
            map_parameters.sparsity_hash_threshold
                = args::get(map_sparsification) * std::numeric_limits<uint64_t>::max();
        }
    } else {
        map_parameters.sparsity_hash_threshold
            = std::numeric_limits<uint64_t>::max();
    }

    if (sketchSize) {
        const int64_t ss = args::get(sketchSize);
        if (ss < 1) {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, sketch size must be at least 1" << std::endl;
            exit(1);
        }
        map_parameters.sketchSize = ss;
    } else {
        map_parameters.sketchSize = 25;
    }

    map_parameters.split = !args::get(no_split);
    map_parameters.mergeMappings = !args::get(no_merge);

    if (segment_length) {
        const int64_t s = mashmap::handy_parameter(args::get(segment_length));

        if (s <= 0) {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, segment length has to be a float value greater than 0." << std::endl;
            exit(1);
        }

        if (s < 100) {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, minimum segment length is required to be >= 100 bp." << std::endl
                      << "[mashmap] This is because Mashmap is not designed for computing short local alignments." << std::endl;
            exit(1);
        }
        map_parameters.segLength = s;
    } else {
        map_parameters.segLength = 5000;
    }

    if (map_pct_identity) {
        if (args::get(map_pct_identity) < 50) {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, minimum nucleotide identity requirement should be >= 50\%." << std::endl;
            exit(1);
        }
        map_parameters.percentageIdentity = (float) (args::get(map_pct_identity)/100.0); // scale to [0,1]
    } else {
        map_parameters.percentageIdentity = skch::fixed::percentage_identity;
    }

    if (map_ani_threshold) {
        map_parameters.ANIDiff = (float) (args::get(map_ani_threshold)/100.0); // scale to [0,1]
    } else {
        map_parameters.ANIDiff = skch::fixed::ANIDiff;
    }

    if (map_ani_threshold_conf) {
        map_parameters.ANIDiffConf = (float) (args::get(map_ani_threshold_conf)/100.0); // scale to [0,1]
    } else {
        map_parameters.ANIDiffConf = skch::fixed::ANIDiffConf;
    }

    map_parameters.stage1_topANI_filter = !args::get(disable_topANI_filter); 
    map_parameters.stage2_full_scan = args::get(stage2_full_scan); 

    if (block_length) {
        const int64_t l = mashmap::handy_parameter(args::get(block_length));

        if (l < 0) {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, min block length has to be a float value greater than or equal to 0." << std::endl;
            exit(1);
        }

        map_parameters.block_length = l;
    } else {
        // Default block length is just the segment size
        map_parameters.block_length = map_parameters.segLength;
    }

    map_parameters.chain = !no_merge;

    if (chain_gap) {
        const int64_t l = mashmap::handy_parameter(args::get(chain_gap));
        if (l < 0) {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, chain gap has to be a float value greater than or equal to 0." << std::endl;
            exit(1);
        }
        map_parameters.chain_gap = l;
    } else {
        // Default chain length is just the segment size / sketch_size
        map_parameters.chain_gap = map_parameters.segLength / map_parameters.sketchSize;
    }

    if (drop_low_map_pct_identity) {
        map_parameters.keep_low_pct_id = false;
    } else {
        map_parameters.keep_low_pct_id = true;
    }

    if (kmer_size) {
        map_parameters.kmerSize = args::get(kmer_size);
    } else {
        // Smaller values of k are more sensitive for divergent genomes, but lose specificity for large
        // genomes due to chance k-mer collisions. However, too large of a k-mer will reduce sensitivity
        // and so choosing the smallest k that avoids chance collisions is recommended.
        /*
        map_parameters.kmerSize = (map_parameters.percentageIdentity >= 0.97 ? 18 :
                                  (map_parameters.percentageIdentity >= 0.9 ? 17 : 15));
        */
        map_parameters.kmerSize = 19;
    }

    if (kmer_pct_threshold) {
        map_parameters.kmer_pct_threshold = args::get(kmer_pct_threshold);
    } else {
        map_parameters.kmer_pct_threshold = 0.001; // in percent! so we keep 99.999% of kmers
    }


    if (thread_count) {
        map_parameters.threads = args::get(thread_count);
    } else {
        map_parameters.threads = 1;
    }


    map_parameters.outFileName = "/dev/stdout";
    yeet_parameters.approx_mapping = true;


    if (num_mappings_for_segments) {
        if (args::get(num_mappings_for_segments) > 0) {
            map_parameters.numMappingsForSegment = args::get(num_mappings_for_segments) ;
        } else {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, the number of mappings to retain for each segment has to be grater than 0." << std::endl;
            exit(1);
        }
    } else {
        map_parameters.numMappingsForSegment = 1;
    }

    if (num_mappings_for_short_seq) {
        if (args::get(num_mappings_for_short_seq) > 0) {
            map_parameters.numMappingsForShortSequence = args::get(num_mappings_for_short_seq);
        } else {
            std::cerr << "[mashmap] ERROR, skch::parseandSave, the number of mappings to retain for each sequence shorter than segment length has to be grater than 0." << std::endl;
            exit(1);
        }
    } else {
        map_parameters.numMappingsForShortSequence = 1;
    }

    //Check if files are valid
    skch::validateInputFiles(map_parameters.querySequences, map_parameters.refSequences);

    temp_file::set_keep_temp(args::get(keep_temp_files));

}

}
