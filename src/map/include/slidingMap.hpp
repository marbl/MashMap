/**
 * @file    slidingMap.hpp
 * @brief   implements ordered map to compute Jaccard
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SLIDING_MAP_HPP 
#define SLIDING_MAP_HPP

#include <vector>
#include <algorithm>
#include <numeric>
#include <map>

//Own includes
#include "map/include/base_types.hpp"

//External includes
#include "assert.hpp"

namespace skch
{
  /**
   * @class     skch::SlideMapper
   * @brief     L1 and L2 mapping stages
   */
  template <typename Q_Info>
    class SlideMapper
    {

      private:

        //Metadata for the minmers saved in sliding ordered map during L2 stage
        struct slidingMapContainerValueType
        {
          hash_t hash;
          strand_t q_strand;
          strand_t strand_vote;
          unsigned int num_before_inc; // Number of hashes between this and the previous (including self)
          bool active;
        };

        //Container type for saving read sketches during L1 and L2 both
        typedef Sketch::MI_Type MinVec_Type;

        typedef Sketch::MIIter_t MIIter_t;

        //reference to query's metadata
        const Q_Info &Q;

        //Define a Not available position marker
        static const offset_t NAPos = std::numeric_limits<offset_t>::max();

        //Ordered map to save unique sketch elements, and associated value as 
        //a pair of its occurrence in the query and the reference
        typedef std::vector<slidingMapContainerValueType> VecType;
        VecType slidingWindowMinhashes;

        //Iterator pointing to the last query minmer that is below rank sketch-size
        typename VecType::iterator pivot;
        typename VecType::size_type pivRank;


      public:

        //Count of shared sketch elements between query and the reference
        //Updated after insert or delete operation on map 
        int sharedSketchElements;

        // Count of strand votes
        int strand_votes;


        //Delete default constructor
        SlideMapper() = delete;

        /**
         * @brief                 constructor
         * @param[in]   Q         query meta data
         */
        SlideMapper(Q_Info &Q_) :
          Q(Q_),
          slidingWindowMinhashes(Q.sketchSize + 1),
          sharedSketchElements(0),
          strand_votes(0)
        {
          this->init();
        }

      private:

        static constexpr bool slidingMapContainerValueType_comp(
            const slidingMapContainerValueType& a, const slidingMapContainerValueType& b) 
        {
          return a.hash < b.hash;
        } 

        /**
         * @brief       Fills map with minimum 's' minmers in the query
         */
        inline void init()
        {
          //Range of sketch in query
          int idx = 0;
          for(auto it = Q.minmerTableQuery.begin(); it != Q.minmerTableQuery.end(); it++)
          {
            this->slidingWindowMinhashes[idx++] = slidingMapContainerValueType {
                it->hash, 
                it->strand, 
                0,  // strand_vote
                1,  // num_before_inc
                0}; // Active (shared)
          }

          // Sort the hashes
          std::sort(slidingWindowMinhashes.begin(), slidingWindowMinhashes.end(), slidingMapContainerValueType_comp);

          //Point pivot to last element in the map
          this->pivot = std::prev(this->slidingWindowMinhashes.end());
          pivRank = slidingWindowMinhashes.size() - 1;
        }

      public:

        inline void insert_minmer(const skch::MinmerInfo& mi)
        {
          // Find where minmer goes in vector
          auto insert_loc = std::lower_bound(
              std::next(slidingWindowMinhashes.begin()), slidingWindowMinhashes.end(),
              mi.hash,
              [](const slidingMapContainerValueType& a, hash_t b) {return a.hash < b;});

          if (insert_loc == slidingWindowMinhashes.end()) 
          {
            return;
          } 

          // If minmer matches, then set curr to active
          if (insert_loc->hash == mi.hash) {
            insert_loc->active = true;
            insert_loc->strand_vote += (insert_loc->q_strand * mi.strand);

            // If minmer matches and is <= pivot, sharedSketchElements++
            if (insert_loc->hash <= pivot->hash) {
              sharedSketchElements++;
              strand_votes += insert_loc->strand_vote; 
            }
          } 
          // Else increment num_before_inc of curr
          else {
            insert_loc->num_before_inc++;
            //  If curr <= pivot then increment pivRank
            if (insert_loc->hash <= pivot->hash) {
              pivRank++;
            }
            //  If pivRank > sketchSize, move piv left and update counters 
            if (pivRank > Q.sketchSize) {
              sharedSketchElements -= pivot->active;
              strand_votes -= pivot->strand_vote; 
              pivRank -= pivot->num_before_inc;
              pivot--;
            }
          }
        }

        /**
         * @brief               delete a minmer from the reference sequence from the map
         * @param[in]   m       reference minmer to remove
         */
        void delete_minmer(const skch::MinmerInfo& mi)
        {
          // Find where minmer goes in vector
          auto insert_loc = std::lower_bound(
              std::next(slidingWindowMinhashes.begin()), slidingWindowMinhashes.end(),
              mi.hash,
              [](const slidingMapContainerValueType& a, hash_t b) {return a.hash < b;});

          if (insert_loc == slidingWindowMinhashes.end()) 
          {
            return;
          } 

          // If minmer matches, then set curr to inactive
          if (insert_loc->hash == mi.hash) {
            // If minmer matches and is <= pivot, sharedSketchElements--
            if (insert_loc->hash <= pivot->hash) {
              sharedSketchElements--;
              strand_votes -= insert_loc->strand_vote;
            }
            insert_loc->active = false;
            insert_loc->strand_vote = 0;
          } 
          // Else decrement num_before_inc of curr
          else {
            insert_loc->num_before_inc--;
            //  If curr <= pivot then decrement pivRank
            if (insert_loc->hash <= pivot->hash) {
              pivRank--;
            }
            //  If pivRank < sketchSize, try to move piv right and update counters 
            if (std::next(pivot) != slidingWindowMinhashes.end()
                && pivRank + std::next(pivot)->num_before_inc <= Q.sketchSize) {
              pivot++;
              sharedSketchElements += pivot->active;
              strand_votes += pivot->strand_vote; 
              pivRank += pivot->num_before_inc;
            }
          }
        }
    };
}

#endif
