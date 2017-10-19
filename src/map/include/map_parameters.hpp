/**
 * @file    map_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_CONFIG_HPP 
#define SKETCH_CONFIG_HPP

#include <vector>

namespace skch
{
  /**
   * @brief   configuration parameters for building sketch
   *          expected to be initialized using command line arguments
   */
  struct Parameters
  {
    int kmerSize;                                     //kmer size for sketching
    int windowSize;                                   //window size used for sketching 
    int segLength;                                    //For split mapping case, this represents the fragment length
                                                      //for noSplit, it represents minimum read length to map
                                                      
    int alphabetSize;                                 //alphabet size
    uint64_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    int filterMode;                                   //filtering mode in mashmap
    int threads;                                      //execution thread count
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
    bool split;                                       //Split read mapping (done if this is true)
  };


  /**
   * @brief     Internal figures not exposed at the command line interface
   */
  namespace fixed
  {
    float pval_cutoff = 1e-03;                        //p-value cutoff for determining window size

    float confidence_interval = 0.75;                 //Confidence interval to relax jaccard cutoff for mapping (0-1)

    float filter_score_best_range = .99;              //mapping score above a certain fraction of best score is 
                                                      //considered good by filtering algorithm

    int max_best_mappings_per_position = 25;          //At a particular position, if algorithm finds more than a certain best 
                                                      //mappings, it doesn't mark them as best anymore
  }
}

#endif
