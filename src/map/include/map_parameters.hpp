/**
 * @file    map_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_CONFIG_HPP
#define SKETCH_CONFIG_HPP

#include <vector>
#include <unordered_set>
#include <filesystem>
namespace stdfs = std::filesystem;

#include "common/ALeS.hpp"

namespace skch
{


struct ales_params {
  uint32_t weight{} ;
  uint32_t seed_count{};
  float similarity{};
  uint32_t region_length{};
};

/**
 * @brief   configuration parameters for building sketch
 *          expected to be initialized using command line arguments
 */
struct Parameters
{
    int kmerSize;                                     //kmer size for sketching
    float kmer_pct_threshold;                         //use only kmers not in the top kmer_pct_threshold %-ile
    int64_t segLength;                                //For split mapping case, this represents the fragment length
                                                      //for noSplit, it represents minimum read length to multimap
    int64_t block_length;                             // minimum (potentially merged) block to keep if we aren't split
    int64_t chain_gap;                                // max distance for 2d range union-find mapping chaining
    int alphabetSize;                                 //alphabet size
    uint64_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    bool stage2_full_scan;                            //Instead of using the best intersection for a given candidate region, compute the minhash for every position in the window
    bool stage1_topANI_filter;                        //Use the ANI filter in stage 1
    float ANIDiff;                                    //ANI distance threshold below best mapping to retain in stage 1 filtering
    float ANIDiffConf;                                //Confidence of stage 1 ANI filtering threshold
    int filterMode;                                   //filtering mode in mashmap
    uint32_t numMappingsForSegment;                   //how many mappings to retain for each segment
    uint32_t numMappingsForShortSequence;             //how many secondary alignments we keep for reads < segLength
    int threads;                                      //execution thread count
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
    stdfs::path saveIndexFilename;                    //output file name of index
    stdfs::path loadIndexFilename;                    //input file name of index
    bool split;                                       //Split read mapping (done if this is true)
    bool skip_self;                                   //skip self mappings
    bool skip_prefix;                                 //skip mappings to sequences with the same prefix
    char prefix_delim;                                //the prefix delimiter
    bool mergeMappings;                               //if we should merge consecutive segment mappings
    bool keep_low_pct_id;                             //true if we should keep mappings whose estimated identity < percentageIdentity
    bool report_ANI_percentage;                       //true if ANI should be in [0,100] as opposed to [0,1] (this is necessary for wfmash


    int sketchSize;
    bool use_spaced_seeds;                            //
    ales_params spaced_seed_params;                   //
    double spaced_seed_sensitivity;                   //
    std::vector<ales::spaced_seed> spaced_seeds;      //
    bool world_minimizers;
    uint64_t sparsity_hash_threshold;                 // keep mappings that hash to <= this value

    bool legacy_output;
    //std::unordered_set<std::string> high_freq_kmers;  //
};


/**
 * @brief     Default values or internal figures not exposed at the command line interface
 */
namespace fixed
{

//float filter_score_best_range = .99;              // mapping score above a certain fraction of best score is
//considered good by filtering algorithm

//int max_best_mappings_per_position = 25;          // At a particular position, if algorithm finds more than a certain best
//mappings, it doesn't mark them as best anymore

double pval_cutoff = 1e-3;                          // p-value cutoff for determining window size
float confidence_interval = 0.95;                   // Confidence interval to relax jaccard cutoff for mapping (0-1)
float percentage_identity = 0.85;                   // Percent identity in the mapping step
float ANIDiff = 0.0;                                // Stage 1 ANI diff threshold
float ANIDiffConf = 0.999;                          // ANI diff confidence
std::string VERSION = "3.0.6";                      // Version of MashMap
}
}

#endif
