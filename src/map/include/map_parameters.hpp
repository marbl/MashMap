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
    int minMatchLength;                               //minimum match length which mashmap should guarantee
    int alphabetSize;                                 //alphabet size
    uint64_t referenceSize;                           //Approximate reference size
    float percentageIdentity;                         //user defined threshold for good similarity
    double p_value;                                   //user defined threshold for p value
    int filterMode;                                   //filtering mode in mashmap
    int threads;                                      //execution thread count
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string outFileName;                          //output file name
    bool split;                                       //Split read mapping (done if this is true)
  };
}

#endif
