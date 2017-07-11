/**
 * @file    filter.hpp
 * @brief   implments the routines to filter mappings
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef FILTER_MAP_HPP 
#define FILTER_MAP_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <zlib.h>  

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"

//External includes

namespace skch
{
  /**
   * @namespace skch::CommonFunc
   * @brief     Implements routines to filter mappings
   */
  namespace Filter
  {

    //Comparators to compare two mappings by their mapping scores
    //Mapping score of a mapping is defined as it's length on the axis times identity
    struct MapScoreCmp
    {
      //Lexographical less than comparison based on mapping length (on query axis) and identity
      static bool less_q(const MappingResult& x, const MappingResult& y) {

        offset_t xlen = x.qlen();
        offset_t ylen = y.qlen();

        return xlen * x.nucIdentity < ylen * y.nucIdentity;
      }

      //comparable (query axis)
      static bool comparable_q(const MappingResult& x, const MappingResult& best) {
        return x.qlen() * x.nucIdentity >= 0.99 * best.qlen() * best.nucIdentity;
      }

      //Lexicographical less than comparison based on mapping length (on reference axis) and identity
      static bool less_r(const MappingResult& x, const MappingResult& y) {

        offset_t xlen = x.rlen();
        offset_t ylen = y.rlen();

        return xlen * x.nucIdentity < ylen * y.nucIdentity;
      }

      //comparable (reference axis)
      static bool comparable_r(const MappingResult& x, const MappingResult& best) {
        return x.rlen() * x.nucIdentity >= 0.99 * best.rlen() * best.nucIdentity;
      }

    };


    /**
     * @namespace skch::filter::query
     * @brief     filter routines (best for query sequence)
     */
    namespace query
    {
      /**
       * @brief                       filter mappings (best for query sequence) 
       * @details                     if any mapping is contained inside another mapping on query axis,
       *                              it is discarded
       * @param[in/out] readMappings  Mappings computed by Mashmap
       */
      template <typename VecIn>
        void filterByContainment(VecIn &readMappings)
        {
          if(readMappings.size() <= 1)
            return;

          //Sort the mappings by query position and longer length
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
                offset_t alen = a.qlen();
                offset_t blen = b.qlen();

                return (a.queryStartPos < b.queryStartPos) || 
                (a.queryStartPos == b.queryStartPos && std::tie(alen, a.nucIdentity) > std::tie(blen, b.nucIdentity));
              });

          //Reference to store longer mapping that may contain smaller ones
          auto longMappingCover = readMappings.begin();

          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            auto thisStartPos = it->queryStartPos;
            auto thisEndPos = it->queryEndPos;

            auto longStartPos = longMappingCover->queryStartPos;
            auto longEndPos = longMappingCover->queryEndPos;

            if(thisEndPos > longEndPos)
            {
              //revise cover mapping
              longMappingCover = it;
            }
            else if(thisEndPos <= longEndPos && !MapScoreCmp::comparable_q(*it, *longMappingCover))
            {
              it->discard = 1;  //contained and lower score
            }
          }

          readMappings.erase( 
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());
        }

     /**
       * @brief                       filter mappings (best for query sequence) 
       * @details                     evaluate best cover mapping for each base pair
       * @param[in/out] readMappings  Mappings computed by Mashmap
       */
      template <typename VecIn>
        void liFilterAlgorithm(VecIn &readMappings)
        {
          if(readMappings.size() <= 1)
            return;

          //Sort the mappings by query position
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
              return (a.queryStartPos < b.queryStartPos);
              });

          //Initially mark all mappings as bad
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ e.discard = 1; });  

          offset_t lastMappingEndPosition;
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e)
              { 
                lastMappingEndPosition = std::max(e.queryEndPos, lastMappingEndPosition); 
              });

          //end position of the best mapping covering the current position of axis being considered
          offset_t bestMappingEndPosition;

          /*
           * subset of mappings that pass the current base pair position
           * each mapping is stored by it's index in the readMappings vector 
           */
          std::list<uint64_t> mappingSubset; 

          auto it = readMappings.begin();

          for(offset_t currentBasePosition = 0; currentBasePosition <= lastMappingEndPosition;)
          {
            //Remove mappings from list that end before currentBasePosition
            mappingSubset.remove_if( [&](uint64_t i){ return readMappings[i].queryEndPos < currentBasePosition; });

            if(mappingSubset.size() == 0)
            {
              //This means there are no mappings that span currentBasePosition
              //Reset bestMappingEndPosition so that next read mappings in the vector are considered
              bestMappingEndPosition = std::numeric_limits<offset_t>::max() - 1;
            }
            else //means mappingSubset.size() > 0
            {
              //Find out best mapping with highest score
              auto mappingIndexWithBestScore = *std::max_element(mappingSubset.begin(), mappingSubset.end(), [&](uint64_t i, uint64_t j){
                  return MapScoreCmp::less_q(readMappings[i], readMappings[j]);
                  });

              auto num_good_mappings = std::count_if(mappingSubset.begin(), mappingSubset.end(), [&](uint64_t i) {
                  return MapScoreCmp::comparable_q(readMappings[i], readMappings[mappingIndexWithBestScore]);
                  });

              if(num_good_mappings <= 25)
              {
                //Unmark the discard flag of good mappings (so that we don't delete them later)
                std::for_each(mappingSubset.begin(), mappingSubset.end(), [&](uint64_t i) { 
                    if( MapScoreCmp::comparable_q(readMappings[i], readMappings[mappingIndexWithBestScore]) )
                    {
                      readMappings[i].discard = -1;
                    }
                  });  
              }
              //else we have too many good mappings over this position, so still ignore these

              //revise bestMappingEndPosition
              bestMappingEndPosition = readMappings[mappingIndexWithBestScore].queryEndPos;
            }

            /*
             *  Advance currentBasePosition ahead
             */

            if (it != readMappings.end())
            {
              auto itStartPos = it->queryStartPos;

              currentBasePosition = std::min(itStartPos, bestMappingEndPosition + 1);

              if( currentBasePosition == itStartPos)
              {
                //Bucket by same query start pos
                auto it_end = std::find_if(it, readMappings.end(), [&](const MappingResult &e){return e.queryStartPos != itStartPos;} ); 

                //Insert mappings in the range [it, it_end) into list

                uint64_t indexOffset = std::distance(readMappings.begin(), it);

                for(uint64_t i = 0; i < std::distance(it, it_end); i++)
                  mappingSubset.push_back(i + indexOffset);

                //advance the iterator
                it = it_end;
              }
            }
            else
            {
              //it now equals readMappings.end() which means we are done pushing all the mappings into list
              currentBasePosition = bestMappingEndPosition + 1;
            }
          }

          readMappings.erase( 
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());
        }

      /**
       * @brief                       filter mappings (best for query sequence) 
       * @param[in/out] readMappings  Mappings computed by Mashmap (post merge step)
       */
      template <typename VecIn>
        void filterMappings(VecIn &readMappings)
        {
          //Remove any mapping contained inside a larger mapping on that axis
          filterByContainment(readMappings);

          //Apply the main filtering algorithm to ensure best mappings across complete axis
          liFilterAlgorithm(readMappings);
        }
    }

    /**
     * @namespace skch::filter::ref
     * @brief     filter routines (best for reference sequence)
     * @details   implementation quite similar to filtering code over query sequence done before
     */
    namespace ref
    {
      //type for saving reference sequence position
      //<seq id, seq position>
      using RefPos_t = std::pair<seqno_t, offset_t>;

      /**
       * @brief                       filter mappings (best for reference sequence) 
       * @details                     if any mapping is contained inside another mapping on reference axis,
       *                              it is discarded
       * @param[in/out] readMappings  Mappings computed by Mashmap
       */
      template <typename VecIn>
        void filterByContainment(VecIn &readMappings)
        {
          if(readMappings.size() <= 1)
            return;

          //Sort the mappings by query position and longer length
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
                offset_t alen = a.rlen();
                offset_t blen = b.rlen();

                return (std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos)) || 
                  (std::tie(a.refSeqId, a.refStartPos) == std::tie(b.refSeqId, b.refStartPos) &&
                   std::tie(alen, a.nucIdentity) > std::tie(blen, b.nucIdentity));
              });

          //Reference to store longer mapping that may contain smaller ones
          auto longMappingCover = readMappings.begin();

          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            RefPos_t thisStartPos = {it->refSeqId, it->refStartPos};
            RefPos_t thisEndPos = {it->refSeqId, it->refEndPos};

            RefPos_t longStartPos = {longMappingCover->refSeqId, longMappingCover->refStartPos};
            RefPos_t longEndPos = {longMappingCover->refSeqId, longMappingCover->refEndPos};

            if(thisEndPos > longEndPos)
            {
              //revise cover mapping
              longMappingCover = it;
            }
            else if(thisEndPos <= longEndPos && !MapScoreCmp::comparable_r(*it, *longMappingCover))
            {
              it->discard = 1;  //contained and lower score
            }
          }

          readMappings.erase( 
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());
        }

      /**
       * @brief                       advance position by one over reference sequence(s)
       * @param[in/out] readMappings  Mappings computed by Mashmap
       * @param[in]     refsketch     reference index class object
       */
      RefPos_t refPosDoPlusOne(const RefPos_t &refPos, const skch::Sketch &refsketch)
      {
        RefPos_t returnPos;

        seqno_t currentSeqId = refPos.first;
        offset_t currentSeqOffSet = refPos.second;

        //if input is max (special reset value in the algorithm), return as it is
        if (currentSeqId == std::numeric_limits<seqno_t>::max() && currentSeqOffSet == std::numeric_limits<offset_t>::max())
          returnPos = refPos;

        //if offset is at the end of reference sequence, shift to next
        else if(currentSeqOffSet == refsketch.metadata[currentSeqId].len - 1)
          returnPos = {currentSeqId + 1, 0};

        else
          returnPos = {currentSeqId, currentSeqOffSet + 1};

        return returnPos;
      }

      /**
       * @brief                       filter mappings (best for reference sequence) 
       * @details                     evaluate best cover mapping for each base pair
       * @param[in/out] readMappings  Mappings computed by Mashmap
       * @param[in]     refsketch     reference index class object
       */
      template <typename VecIn>
        void liFilterAlgorithm(VecIn &readMappings, const skch::Sketch &refsketch)
        {
          if(readMappings.size() <= 1)
            return;

          //Sort the mappings by query position
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
                return std::tie(a.refSeqId, a.refStartPos) < std::tie(b.refSeqId, b.refStartPos);
              });

          //Initially mark all mappings as bad
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ e.discard = 1; });  

          RefPos_t lastMappingEndPosition;
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e)
              { 
                lastMappingEndPosition = std::max(std::make_pair(e.refSeqId, e.refEndPos), lastMappingEndPosition); 
              });

          //end position of the best mapping covering the current position of axis being considered
          RefPos_t bestMappingEndPosition;

          /*
           * subset of mappings that pass the current base pair position
           * each mapping is stored by it's index in the readMappings vector 
           */
          std::list<uint64_t> mappingSubset; 

          auto it = readMappings.begin();

          for(RefPos_t currentBasePosition = {0,0}; currentBasePosition <= lastMappingEndPosition;)
          {
            //Remove mappings from list that end before currentBasePosition
            mappingSubset.remove_if( [&](uint64_t i){ return std::make_pair(readMappings[i].refSeqId, readMappings[i].refEndPos) < currentBasePosition; });

            if (mappingSubset.size() == 0)
            {
              //This means there are no mappings that span currentBasePosition
              //Reset bestMappingEndPosition so that next read mappings in the vector are considered
              bestMappingEndPosition = {std::numeric_limits<seqno_t>::max(), std::numeric_limits<offset_t>::max()};
            }
            else //means mappingSubset.size() > 0
            {
              //Find out best mapping with highest score
              auto mappingIndexWithBestScore = *std::max_element(mappingSubset.begin(), mappingSubset.end(), [&](uint64_t i, uint64_t j){
                  return MapScoreCmp::less_r(readMappings[i], readMappings[j]);
                  });

              auto num_good_mappings = std::count_if(mappingSubset.begin(), mappingSubset.end(), [&](uint64_t i) {
                  return MapScoreCmp::comparable_r(readMappings[i], readMappings[mappingIndexWithBestScore]);
                  });

              if(num_good_mappings <= 25)
              {
                //Unmark the discard flag of good mappings (so that we don't delete them later)
                std::for_each(mappingSubset.begin(), mappingSubset.end(), [&](uint64_t i) { 
                    if( MapScoreCmp::comparable_r(readMappings[i], readMappings[mappingIndexWithBestScore]) )
                    {
                      readMappings[i].discard = -1;
                    }
                  });  
              }
              //else we have too many good mappings over this position, so still ignore these

              //revise bestMappingEndPosition
              bestMappingEndPosition = {readMappings[mappingIndexWithBestScore].refSeqId, readMappings[mappingIndexWithBestScore].refEndPos};
            }

            /*
             *  Advance currentBasePosition ahead
             */

            if (it != readMappings.end())
            {
              RefPos_t itStartPos = {it->refSeqId, it->refStartPos};

              currentBasePosition = std::min( itStartPos, refPosDoPlusOne(bestMappingEndPosition, refsketch) );

              if( currentBasePosition == itStartPos)
              {
                //Bucket by same reference start pos
                auto it_end = std::find_if(it, readMappings.end(), [&](const MappingResult &e){return std::make_pair(e.refSeqId, e.refStartPos) != itStartPos;} ); 

                //Insert mappings in the range [it, it_end) into list

                uint64_t indexOffset = std::distance(readMappings.begin(), it);

                for(uint64_t i = 0; i < std::distance(it, it_end); i++)
                  mappingSubset.push_back(i + indexOffset);

                //advance the iterator
                it = it_end;
              }
            }
            else
            {
              //it now equals readMappings.end() which means we are done pushing all the mappings into list
              currentBasePosition = refPosDoPlusOne(bestMappingEndPosition, refsketch);
            }
          }

          readMappings.erase( 
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());
        }

      /**
       * @brief                       filter mappings (best for query sequence) 
       * @param[in/out] readMappings  Mappings computed by Mashmap (post merge step)
       * @param[in]     refsketch     reference index class object, used to determine ref sequence lengths
       */
      template <typename VecIn>
        void filterMappings(VecIn &readMappings, const skch::Sketch &refsketch)
        {
          //Remove any mapping contained inside a larger mapping on that axis
          filterByContainment(readMappings);

          //Apply the main filtering algorithm to ensure best mappings across complete axis
          liFilterAlgorithm(readMappings, refsketch);
        }
    }

  }
}

#endif
