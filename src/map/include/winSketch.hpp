/**
 * @file    winSketch.hpp
 * @brief   routines to index the reference 
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WIN_SKETCH_HPP 
#define WIN_SKETCH_HPP

#include <algorithm>
#include <cassert>
#include <map>
#include <unordered_map>
#include <unordered_set>
#include <vector>
//#include <zlib.h>

//Own includes
#include "csv.h"
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/ThreadPool.hpp"

//External includes
#include "common/murmur3.h"
#include "common/prettyprint.hpp"

//#include "common/sparsehash/dense_hash_map"
//#include "common/parallel-hashmap/parallel_hashmap/phmap.h"
//#include <abseil-cpp/absl/container/flat_hash_map.h>
//#include <common/sparse-map/include/tsl/sparse_map.h>
#include <common/robin-hood-hashing/robin_hood.h>

#include "common/seqiter.hpp"

#include "assert.hpp"

namespace skch
{
  /**
   * @class     skch::Sketch
   * @brief     sketches and indexes the reference (subject sequence)
   * @details  
   *            1.  Minmers are computed in streaming fashion
   *                Computing minmers is using double ended queue which gives
   *                O(reference size) complexity
   *                Algorithm described here:
   *                https://people.cs.uct.ac.za/~ksmith/articles/sliding_window_minimum.html
   *
   *            2.  Index hashes into appropriate format to enable fast search at L1 mapping stage
   */
  class Sketch
    {
      //private members
    
      //algorithm parameters
      const skch::Parameters &param;

      //Minmers that occur this or more times will be ignored (computed based on percentageThreshold)
      int freqThreshold = std::numeric_limits<int>::max();

      //Set of frequent seeds to be ignored
      std::unordered_set<hash_t> frequentSeeds;

      //Make the default constructor private, non-accessible
      Sketch();

      public:

      typedef std::vector< MinmerInfo > MI_Type;
      using MIIter_t = MI_Type::const_iterator;

      //Keep sequence length, name that appear in the sequence (for printing the mappings later)
      std::vector< ContigInfo > metadata;

      /*
       * Keep the information of what sequences come from what file#
       * Example [a, b, c] implies 
       *  file 0 contains 0 .. a-1 sequences
       *  file 1 contains a .. b-1 
       *  file 2 contains b .. c-1
       */
      std::vector< int > sequencesByFileInfo;

      //Index for fast seed lookup (unordered_map)
      /*
       * [minmer #1] -> [pos1, pos2, pos3 ...]
       * [minmer #2] -> [pos1, pos2...]
       * ...
       */
      //using MI_Map_t = google::dense_hash_map< MinmerMapKeyType, MinmerMapValueType >;
      //using MI_Map_t = phmap::flat_hash_map< MinmerMapKeyType, MinmerMapValueType >;
      //using MI_Map_t = absl::flat_hash_map< MinmerMapKeyType, MinmerMapValueType >;
      //using MI_Map_t = tsl::sparse_map< MinmerMapKeyType, MinmerMapValueType >;
      using MI_Map_t = robin_hood::unordered_flat_map< MinmerMapKeyType, MinmerMapValueType >;
      MI_Map_t minmerPosLookupIndex;
      MI_Type minmerIndex;

      private:

      /**
       * Keep list of minmers, sequence# , their position within seq , here while parsing sequence 
       * Note : position is local within each contig
       * Hashes saved here are non-unique, ordered as they appear in the reference
       */

      //Frequency histogram of minmers
      //[... ,x -> y, ...] implies y number of minmers occur x times
      std::map<int, int> minmerFreqHistogram;

      public:

      /**
       * @brief   constructor
       *          also builds, indexes the minmer table
       */
      Sketch(const skch::Parameters &p) 
        :
          param(p) {
            this->build();
            if (param.loadIndexFilename == "") {
              if (param.saveIndexFilename != "") {
                this->saveIndex();
              }
            }
            this->index();
            this->computeFreqHist();
            this->computeFreqSeedSet();
            this->dropFreqSeedSet();
          }

      private:

      /**
       * @brief     build the sketch table
       * @details   compute and save minmers from the reference sequence(s)
       *            assuming a fixed window size
       */
      void build()
      {
        //sequence counter while parsing file
        seqno_t seqCounter = 0;

        //Create the thread pool 
        ThreadPool<InputSeqContainer, MI_Type> threadPool( [this](InputSeqContainer* e) {return buildHelper(e);}, param.threads);
        if (param.loadIndexFilename != "") {
          this->loadIndex();
        }

        for(const auto &fileName : param.refSequences)
        {

#ifdef DEBUG
        std::cerr << "[mashmap::skch::Sketch::build] building minmer index for " << fileName << std::endl;
#endif

        seqiter::for_each_seq_in_file(
            fileName,
            [&](const std::string& seq_name,
                const std::string& seq) {
                // todo: offset_t is an 32-bit integer, which could cause problems
                offset_t len = seq.length();

                //Save the sequence name
                metadata.push_back( ContigInfo{seq_name, len} );

                //Is the sequence too short?
                if(len < param.kmerSize)
                {
#ifdef DEBUG
                    std::cerr << "WARNING, skch::Sketch::build, found an unusually short sequence relative to kmer" << std::endl;
#endif
                }
                else
                {
                  if (param.loadIndexFilename == "") {
                    threadPool.runWhenThreadAvailable(new InputSeqContainer(seq, seq_name, seqCounter));
                    
                    //Collect output if available
                    while ( threadPool.outputAvailable() )
                        this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
                  }
                }
                seqCounter++;
            });

          sequencesByFileInfo.push_back(seqCounter);
        }

        if (param.loadIndexFilename == "") {
          //Collect remaining output objects
          while ( threadPool.running() )
            this->buildHandleThreadOutput(threadPool.popOutputWhenAvailable());
        }

        std::cerr << "[mashmap::skch::Sketch::build] minmers picked from reference = " << minmerIndex.size() << std::endl;

      }

      /**
       * @brief               function to compute minmers given input sequence object
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @return              output object containing the mappings
       */
      MI_Type* buildHelper(InputSeqContainer *input)
      {
        MI_Type* thread_output = new MI_Type();

        //Compute minmers in reference sequence
        skch::CommonFunc::addMinmers(
                *thread_output, 
                &(input->seq[0u]), 
                input->len, 
                param.kmerSize, 
                param.segLength, 
                param.alphabetSize, 
                param.sketchSize,
                input->seqCounter);

        return thread_output;
      }

      /**
       * @brief                 routine to handle thread's local minmer index
       * @param[in] output      thread local minmer output
       */
      void buildHandleThreadOutput(MI_Type* output)
      {
        this->minmerIndex.insert(this->minmerIndex.end(), output->begin(), output->end());
        delete output;
      }


      /**
       * @brief  Save index for quick loading
       */
      void saveIndex() 
      {
        std::ofstream outStream;
        outStream.open(param.saveIndexFilename, std::fstream::out);
        outStream << "seqId" << "\t" << "strand" << "\t" << "start" << "\t" << "end" << "\t" << "hash" << std::endl;
        for (auto& mi : this->minmerIndex) {
          outStream << mi.seqId << "\t" << std::to_string(mi.strand) << "\t" << mi.wpos << "\t" << mi.wpos_end << "\t" << mi.hash << std::endl;
        }
        outStream.close(); 
      }

      /**
       * @brief Load index from file
       */
      void loadIndex() 
      {
        io::CSVReader<5, io::trim_chars<' '>, io::no_quote_escape<'\t'>> inReader(param.loadIndexFilename);
        inReader.read_header(io::ignore_missing_column, "seqId", "strand", "start", "end", "hash");
        hash_t hash;
        offset_t start, end;
        strand_t strand;
        seqno_t seqId;
        while (inReader.read_row(seqId, strand, start, end, hash))
        {
          this->minmerIndex.push_back(MinmerInfo {hash, seqId, start, end, strand});
        }
      }

      /**
       * @brief   build the index for fast lookups using minmer table
       */
      void index()
      {
        //Parse all the minmers and push into the map
        //minmerPosLookupIndex.set_empty_key(0);

        for(auto &e : minmerIndex)
        {
          // [hash value -> info about minmer]
          minmerPosLookupIndex[e.hash].push_back(e);

        }

        std::cerr << "[mashmap::skch::Sketch::index] unique minmers = " << minmerPosLookupIndex.size() << std::endl;
      }

      /**
       * @brief   report the frequency histogram of minmers using position lookup index
       *          and compute which high frequency minmers to ignore
       */
      void computeFreqHist()
      {
          if (!minmerPosLookupIndex.empty()) {
              //1. Compute histogram

              for (auto &e : this->minmerPosLookupIndex)
                  this->minmerFreqHistogram[e.second.size()] += 1;

              std::cerr << "[mashmap::skch::Sketch::computeFreqHist] Frequency histogram of minmers = "
                        << *this->minmerFreqHistogram.begin() << " ... " << *this->minmerFreqHistogram.rbegin()
                        << std::endl;

              //2. Compute frequency threshold to ignore most frequent minmers

              int64_t totalUniqueMinmers = this->minmerPosLookupIndex.size();
              int64_t minmerToIgnore = totalUniqueMinmers * param.kmer_pct_threshold / 100;

              int64_t sum = 0;

              //Iterate from highest frequent minmers
              for (auto it = this->minmerFreqHistogram.rbegin(); it != this->minmerFreqHistogram.rend(); it++) {
                  sum += it->second; //add frequency
                  if (sum < minmerToIgnore) {
                      this->freqThreshold = it->first;
                      //continue
                  } else if (sum == minmerToIgnore) {
                      this->freqThreshold = it->first;
                      break;
                  } else {
                      break;
                  }
              }

              if (this->freqThreshold != std::numeric_limits<int>::max())
                  std::cerr << "[mashmap::skch::Sketch::computeFreqHist] With threshold " << this->param.kmer_pct_threshold
                            << "\%, ignore minmers occurring >= " << this->freqThreshold << " times during lookup."
                            << std::endl;
              else
                  std::cerr << "[mashmap::skch::Sketch::computeFreqHist] With threshold " << this->param.kmer_pct_threshold
                            << "\%, consider all minmers during lookup." << std::endl;
          } else {
              std::cerr << "[mashmap::skch::Sketch::computeFreqHist] No minmers." << std::endl;
          }
      }

      public:

      /**
       * @brief               search hash associated with given position inside the index
       * @details             if MIIter_t iter is returned, than *iter's wpos >= winpos
       * @param[in]   seqId
       * @param[in]   winpos
       * @return              iterator to the minmer in the index
       */

      /**
       * @brief                 check if iterator points to index end
       * @param[in]   iterator
       * @return                boolean value
       */
      bool isMinmerIndexEnd(const MIIter_t &it) const
      {
        return it == this->minmerIndex.end();
      }

      /**
       * @brief     Return end iterator on minmerIndex
       */
      MIIter_t getMinmerIndexEnd() const
      {
        return this->minmerIndex.end();
      }

      int getFreqThreshold() const
      {
        return this->freqThreshold;
      }

      void computeFreqSeedSet()
      {
        for(auto &e : this->minmerPosLookupIndex) {
          if (e.second.size() >= this->freqThreshold) {
            this->frequentSeeds.insert(e.first);
          }
        }
      }

      void dropFreqSeedSet()
      {
        this->minmerIndex.erase(
          std::remove_if(minmerIndex.begin(), minmerIndex.end(), [&] 
            (auto& mi) {return this->frequentSeeds.find(mi.hash) != this->frequentSeeds.end();}
          ), minmerIndex.end()
        );
      }

      bool isFreqSeed(hash_t h) const
      {
        return frequentSeeds.find(h) != frequentSeeds.end();
      }

    }; //End of class Sketch
} //End of namespace skch

#endif
