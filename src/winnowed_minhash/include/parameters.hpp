/**
 * @file    winmin_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WINMINHASH_CONFIG_HPP 
#define WINMINHASH_CONFIG_HPP

#include <vector>

namespace winmin
{
  /**
   * @brief   configuration parameters for building sketch
   *          expected to be initialized using command line arguments
   */
  struct Parameters
  {
    int kmerSize;                                     //kmer size for sketching
    int windowSize;                                   //window size used for sketching 
    std::string refSequence;                          //reference sequence
    std::string querySequence;                        //query sequence
  };
}

#endif
