/**
 * @file    filter.hpp
 * @brief   implments the routines to filter mappings
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef FILTER_MAP_HPP 
#define FILTER_MAP_HPP

#include <vector>
#include <algorithm>
#include <set>
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
    /**
     * @namespace skch::filter::query
     * @brief     filter routines (best for query sequence)
     */
    namespace query
    {
      //helper functions for executing plane sweep over query sequence
      struct Helper
      {
        MappingResultsVector_t &vec;

        Helper(MappingResultsVector_t &v) : vec(v) {}

        //Greater than comparison by score and begin position
        //used to define order in BST
        bool operator ()(const int &x, const int &y) const {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = vec[x].qlen() * vec[x].nucIdentity;
          auto y_score = vec[y].qlen() * vec[y].nucIdentity;

          return std::tie(x_score, vec[x].queryStartPos) > std::tie(y_score, vec[y].queryStartPos);
        }

        //Greater than comparison by score
        bool greater_score(const int &x, const int &y) {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = vec[x].qlen() * vec[x].nucIdentity;
          auto y_score = vec[y].qlen() * vec[y].nucIdentity;

          return x_score > y_score;
        }

        /*
         * @brief                         mark the mappings with maximum score as good (on query seq)
         * @tparam          Type          std::set type to save mappings (sweep line status container)
         * @param[in/out]   L             container with mappings
         */
        template <typename Type>
          inline void markGood(Type &L)
          {
            //first segment in the set order
            auto beg = L.begin();

            for(auto it = L.begin(); it != L.end(); it++)
            {
              if(this->greater_score(*beg, *it) || vec[*it].discard == 0) 
                break;

              vec[*it].discard = 0;
            }
          }
      };

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

          //Initially mark all mappings as bad
          //Maintain the order of this vector till end of this function
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ e.discard = 1; });  

          //Initialize object of Helper struct
          Helper obj (readMappings);

          //Plane sweep status
          //binary search tree of segment ids, ordered by their scores
          std::set <int, Helper> bst (obj);

          //Event point schedule
          //vector of triplets <position, event type, segment id>
          typedef std::tuple<offset_t, int, int> eventRecord_t;
          std::vector <eventRecord_t>  eventSchedule (2*readMappings.size());  

          for(int i = 0; i < readMappings.size(); i++)
          {
            eventSchedule.emplace_back (readMappings[i].queryStartPos, event::BEGIN, i);
            eventSchedule.emplace_back (readMappings[i].queryEndPos + 1, event::END, i);
          }

          std::sort(eventSchedule.begin(), eventSchedule.end());

          //Execute the plane sweep algorithm
          for(auto it = eventSchedule.begin(); it!= eventSchedule.end();)
          {
            //Find events that correspond to current position
            auto it2 = std::find_if(it, eventSchedule.end(), [&](const eventRecord_t &e)
                                    {
                                      return std::get<0>(e) != std::get<0>(*it);
                                    });

            //update sweep line status by adding/removing segments
            std::for_each(it, it2, [&](const eventRecord_t &e)
                                    {
                                      if (std::get<1>(e) == event::BEGIN) 
                                        bst.insert (std::get<2>(e));
                                      else 
                                        bst.erase (std::get<2>(e));
                                    });

            //mark mappings as good
            obj.markGood(bst);

            it = it2;
          }

          //Remove bad mappings
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
          //Apply the main filtering algorithm to ensure best mappings across complete axis
          liFilterAlgorithm(readMappings);
        }
    } //End of query namespace

    namespace ref
    {

      //Event point schedule
      //vector of triplets <ref sequence id, ref seq offset, event type, segment id>
      typedef std::tuple <seqno_t, offset_t, int, int> eventRecord_t;

      //helper functions for executing plane sweep over reference sequence
      struct Helper
      {
        MappingResultsVector_t &vec;

        Helper(MappingResultsVector_t &v) : vec(v) {}

        //Greater than comparison by score and begin position
        //used to define order in BST
        bool operator ()(const int &x, const int &y) const {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = vec[x].rlen() * vec[x].nucIdentity;
          auto y_score = vec[y].rlen() * vec[y].nucIdentity;

          return std::tie(x_score, vec[x].refStartPos) > std::tie(y_score, vec[y].refStartPos);
        }

        //Greater than comparison by score
        bool greater_score(const int &x, const int &y) {

          assert(x < vec.size());
          assert(y < vec.size());

          auto x_score = vec[x].rlen() * vec[x].nucIdentity;
          auto y_score = vec[y].rlen() * vec[y].nucIdentity;

          return x_score > y_score;
        }

        /**
         * @brief                         mark the mappings with maximum score as good (on query seq)
         * @tparam          Type          std::set type to save mappings (sweep line status container)
         * @param[in/out]   L             container with mappings
         */
        template <typename Type>
          inline void markGood(Type &L)
          {
            //first segment in the set order
            auto beg = L.begin();

            for(auto it = L.begin(); it != L.end(); it++)
            {
              if(this->greater_score(*beg, *it) || vec[*it].discard == 0) 
                break;

              vec[*it].discard = 0;
            }
          }

        /**
         * @brief                       advance position by one over reference sequence(s)
         * @param[in/out] eventRecord   event record containing end point of segment
         * @param[in]     refsketch     reference index class object
         */
        void refPosDoPlusOne(eventRecord_t &eventRecord, const skch::Sketch &refsketch)
        {
          seqno_t currentSeqId = std::get<0>(eventRecord);
          offset_t currentSeqOffSet = std::get<1>(eventRecord);

          //if offset is at the end of reference sequence, shift to next
          if(currentSeqOffSet == refsketch.metadata[currentSeqId].len - 1)
          {
            std::get<0>(eventRecord) += 1;    //shift id by 1
            std::get<1>(eventRecord) = 0;
          }
          else
            std::get<1>(eventRecord) += 1;    //shift offset by 1:
        }
      };

      /**
       * @brief                       filter mappings (best for reference sequence) 
       * @param[in/out] readMappings  Mappings computed by Mashmap (post merge step)
       * @param[in]     refsketch     reference index class object, used to determine ref sequence lengths
       */
      template <typename VecIn>
        void filterMappings(VecIn &readMappings, const skch::Sketch &refsketch)
        {
          if(readMappings.size() <= 1)
            return;

          //Initially mark all mappings as bad
          //Maintain the order of this vector till end of this function
          std::for_each(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ e.discard = 1; });  

          //Initialize object of Helper struct
          Helper obj (readMappings);

          //Plane sweep status
          //binary search tree of segment ids, ordered by their scores
          std::set <int, Helper> bst (obj);

          //Event point schedule
          //vector of triplets <position, event type, segment id>
          std::vector <eventRecord_t>  eventSchedule (2*readMappings.size());  

          for(int i = 0; i < readMappings.size(); i++)
          {
            eventSchedule.emplace_back (readMappings[i].refSeqId, readMappings[i].refStartPos, event::BEGIN, i);

            eventRecord_t endEvent = std::make_tuple(readMappings[i].refSeqId, readMappings[i].refEndPos, event::END, i);
            //add one to above coordinate
            obj.refPosDoPlusOne(endEvent, refsketch);  
            eventSchedule.push_back (endEvent);
          }

          std::sort(eventSchedule.begin(), eventSchedule.end());

          //Execute the plane sweep algorithm
          for(auto it = eventSchedule.begin(); it!= eventSchedule.end();)
          {
            //Find events that correspond to current position
            auto it2 = std::find_if(it, eventSchedule.end(), [&](const eventRecord_t &e)
                                    {
                                      return std::tie(std::get<0>(e), std::get<1>(e)) != std::tie(std::get<0>(*it), std::get<1>(*it));
                                    });

            //update sweep line status by adding/removing segments
            std::for_each(it, it2, [&](const eventRecord_t &e)
                                    {
                                      if (std::get<2>(e) == event::BEGIN) 
                                        bst.insert (std::get<3>(e));
                                      else 
                                        bst.erase (std::get<3>(e));
                                    });

            //mark mappings as good
            obj.markGood(bst);

            it = it2;
          }

          //Remove bad mappings
          readMappings.erase( 
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());
        } 
    } //End of reference namespace
  }
}

#endif
