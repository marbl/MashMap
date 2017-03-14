/**
 * @file    computeJaccard.hpp
 * @brief   compute jaccard using winnowed-minhash estimator
 *          between 2 given sequences
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WINMIN_JACCARD_HPP 
#define WINMIN_JACCARD_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <zlib.h>  

//Own includes
#include "winnowed_minhash/include/parameters.hpp"
#include "base_types.hpp"
#include "winSketch.hpp"
#include "commonFunc.hpp"
#include "map_stats.hpp"

//External includes

namespace winmin
{
  /**
   * @class     winmin::WinMin
   * @brief     L1 and L2 mapping stages
   */
  class WinMin
  {
    private:

      //algorithm parameters
      const winmin::Parameters &param;

      //Container type for saving read sketch elements
      typedef skch::Sketch::MI_Type MinVec_Type;

      typedef skch::Sketch::MIIter_t MIIter_t;

      //Borrow type definitions from map
      typedef skch::offset_t offset_t;
      typedef skch::seqno_t seqno_t;

    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] refSketch   reference sketch
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      WinMin(const winmin::Parameters &p) :
        param(p)
      {
        this->computeWinminJaccard();
      }

    private:

      /**
       * @brief   compute Jaccard similarity using Win Minhash estimator
       */
      void computeWinminJaccard()
      {
        std::string refSequence = param.refSequence;
        std::string querySequence = param.querySequence;

        MinVec_Type minimizerTableQuery_r, minimizerTableQuery_q;
        offset_t len_r, len_q;

        this->computeMinimizers(param.refSequence, minimizerTableQuery_r, len_r);
        this->computeMinimizers(param.querySequence, minimizerTableQuery_q, len_q);

        std::cout << "Length of the first sequence from reference file = " << len_r << "\n";
        std::cout << "Length of the first sequence from query file = " << len_q << "\n";

        this->jaccardFromMinimizers(minimizerTableQuery_r, minimizerTableQuery_q);
      }

      /**
       * @brief   compute minimizers from a sequence in filename
       */
      void computeMinimizers(const std::string &filename, 
          MinVec_Type &vec, offset_t &len)
      {
        FILE *file = fopen(filename.c_str(), "r");
        gzFile fp = gzdopen(fileno(file), "r");
        kseq_t *seq = kseq_init(fp);
        len = kseq_read(seq);

        skch::CommonFunc::addMinimizers(vec, seq, param.kmerSize, param.windowSize, 4);

        //Close the input file
        kseq_destroy(seq);  
        gzclose(fp);  
      }

      /**
       * @brief   Estimate jaccard similarity using minimizers
       */
      void jaccardFromMinimizers(MinVec_Type &minimizerTableQuery_r, MinVec_Type &minimizerTableQuery_q)
      {
        //Sketch size
        int s;

        std::sort(minimizerTableQuery_q.begin(), minimizerTableQuery_q.end(), skch::MinimizerInfo::lessByHash);
        auto uniqEndIter = std::unique(minimizerTableQuery_q.begin(), minimizerTableQuery_q.end(), skch::MinimizerInfo::equalityByHash);

        s = std::distance(minimizerTableQuery_q.begin(), uniqEndIter);

        /*
         * [Hash value] -> [Existence : 0 (Query), 1 (Ref), 2 (both)]
         */
        typedef std::map< skch::hash_t, int > MapType;
        MapType map;

        int shared_sketch_elements, total_sketch_elements;

        for(auto &e : minimizerTableQuery_q)
          map[e.hash] = 0;                    //Insert as hash from query

        for(auto &e : minimizerTableQuery_r)
          if(map.find(e.hash) == map.end())
            map[e.hash] = 1;                  //Insert as hash from reference           
          else if(map[e.hash] == 0)
            map[e.hash] = 2;                  //Update as hash from both

          //Iterate over map
          int uniqueHashes = 0, sharedHashes = 0;

          for (auto it = map.cbegin(); it != map.cend(); ++it)
          {
            uniqueHashes += 1;

            if(it->second == 2)
              sharedHashes += 1;

            if(uniqueHashes == s)
              break;
          }

          std::cout << "Sketch size = " << s << "\n";
          std::cout << "Shared sketch elements = " << sharedHashes << "\n";
          std::cout << "Jaccard = " << sharedHashes * 1.0/s << "\n";
          std::cout << "Estimated identity = " << 100.0 - 100.0 * skch::Stat::j2md(sharedHashes * 1.0/s, param.kmerSize) << "\n";
      }
  };

}

#endif
