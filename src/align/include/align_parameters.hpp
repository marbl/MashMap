/**
 * @file    align_parameters.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef ALIGN_PARAMETERS_HPP 
#define ALIGN_PARAMETERS_HPP

#include <vector>

namespace align
{
  /**
   * @brief   parameters for generating mashmap alignments
   */
  struct Parameters
  {
    int threads;                                      //execution thread count
    float percentageIdentity;                         //user defined threshold for good similarity
    int bandwidth;                                    //bandwidth cap for edlib
    std::vector<std::string> refSequences;            //reference sequence(s)
    std::vector<std::string> querySequences;          //query sequence(s)
    std::string mashmapPafFile;                       //mashmap paf mapping file
    std::string samOutputFile;                        //sam output file name
  };
}

#endif
