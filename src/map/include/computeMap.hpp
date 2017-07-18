/**
 * @file    computeMap.hpp
 * @brief   implments the sequence mapping logic
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef SKETCH_MAP_HPP 
#define SKETCH_MAP_HPP

#include <vector>
#include <algorithm>
#include <unordered_map>
#include <fstream>
#include <zlib.h>  
#include <cassert>

//Own includes
#include "map/include/base_types.hpp"
#include "map/include/map_parameters.hpp"
#include "map/include/commonFunc.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/slidingMap.hpp"
#include "map/include/MIIteratorL2.hpp"
#include "map/include/filter.hpp"

//External includes

namespace skch
{
  /**
   * @class     skch::Map
   * @brief     L1 and L2 mapping stages
   */
  class Map
  {
    public:

      //Type for Stage L1's predicted candidate location
      struct L1_candidateLocus_t 
      {
        seqno_t seqId;                    //sequence id where read is mapped

        /* read could be mapped with its begin location
         * from [rangeStartPos, rangeEndPos]
         */
        offset_t rangeStartPos;
        offset_t rangeEndPos;  
      };

      //Type for Stage L2's predicted mapping coordinate within each L1 candidate
      struct L2_mapLocus_t 
      {
        seqno_t seqId;                    //sequence id where read is mapped
        offset_t meanOptimalPos;          //Among multiple consecutive optimal positions, save the avg.
        Sketch::MIIter_t optimalStart;    //optimal start mapping position (begin iterator)
        Sketch::MIIter_t optimalEnd;      //optimal end mapping position (end iterator) 
        int sharedSketchSize;             //count of shared sketch elements
      };

    private:

      //algorithm parameters
      const skch::Parameters &param;

      //reference sketch
      const skch::Sketch &refSketch;

      //Container type for saving read sketches during L1 and L2 both
      typedef Sketch::MI_Type MinVec_Type;

      typedef Sketch::MIIter_t MIIter_t;

      //Custom function for post processing the results, by default does nothing 
      typedef std::function< void(const MappingResult&) > PostProcessResultsFn_t;
      PostProcessResultsFn_t processMappingResults;

      //Container to store query sequence name and length
      //used only if one-to-one filtering is ON
      std::vector<ContigInfo> qmetadata; 

      //Count of reads to parse before processing in parallel
      int bufferSize;

    public:

      /**
       * @brief                 constructor
       * @param[in] p           algorithm parameters
       * @param[in] refSketch   reference sketch
       * @param[in] f           optional user defined custom function to post process the reported mapping results
       */
      Map(const skch::Parameters &p, const skch::Sketch &refsketch,
          PostProcessResultsFn_t f = nullptr) :
        param(p),
        refSketch(refsketch),
        processMappingResults(f)
    {
      //Keep the size of input buffer 128 times thread count
      bufferSize = 128 * param.threads;

      this->mapQuery();
    }

    private:

      /**
       * @brief   parse over sequences in query file and map each on the reference
       */
      void mapQuery()
      {
        //Count of reads mapped by us
        //Some reads are dropped because of short length
        seqno_t totalReadsPickedForMapping = 0;
        seqno_t totalReadsMapped = 0;
        seqno_t seqCounter = 0;     

        std::ofstream outstrm(param.outFileName);
        MappingResultsVector_t allReadMappings;  //Aggregate mapping results for the complete run


        std::vector<MapModuleInput> inputBuffer(this->bufferSize);
        std::vector<MapModuleOutput> outputBuffer(this->bufferSize);
        int pendingWorkItems = 0;

        for(const auto &fileName : param.querySequences)
        {
          //Open the file using kseq
          FILE *file = fopen(fileName.c_str(), "r");
          gzFile fp = gzdopen(fileno(file), "r");
          kseq_t *seq = kseq_init(fp);

#ifdef DEBUG
          std::cout << "INFO, skch::Map::mapQuery, mapping reads in " << fileName << std::endl;
#endif

          //size of sequence
          offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {

            if (param.filterMode == filter::ONETOONE)
              qmetadata.push_back( ContigInfo{seq->name.s, (offset_t) seq->seq.l} );

            //Is the read too short?
            if(len < param.windowSize || len < param.kmerSize || len < param.minMatchLength)
            {

#ifdef DEBUG
              std::cout << "WARNING, skch::Map::mapQuery, read is not long enough for mapping" << std::endl;
#endif

              seqCounter++;
              continue;
            }
            else 
            {
              totalReadsPickedForMapping++;

              inputBuffer[pendingWorkItems].init(seq->seq.s, seq->name.s, len, seqCounter);
              pendingWorkItems++;

              //If buffer size is full, process the work
              if(pendingWorkItems == this->bufferSize)
              {
                this->processMappingWorkInParallel(pendingWorkItems, inputBuffer, outputBuffer,
                    allReadMappings, totalReadsMapped, outstrm);

                //Rest pendingWorkItems counter
                pendingWorkItems = 0;
              }
            }

            seqCounter++;
          } //Finish reading query input file

          //Close the input file
          kseq_destroy(seq);  
          gzclose(fp);  
        }

        this->processMappingWorkInParallel(pendingWorkItems, inputBuffer, outputBuffer,
            allReadMappings, totalReadsMapped, outstrm);

        //Filter over reference axis and report the mappings
        if (param.filterMode == filter::ONETOONE)
        {
          skch::Filter::ref::filterMappings(allReadMappings, this->refSketch);
          reportReadMappings(allReadMappings, outstrm);
        }

        std::cout << "INFO, skch::Map::mapQuery, [count of mapped reads, reads qualified for mapping, total input reads] = [" << totalReadsMapped << ", " << totalReadsPickedForMapping << ", " << seqCounter << "]" << std::endl;

      }

      /**
       * @brief                           schedule mapping tasks and output handling inside openmp runtime system
       * @param[in]   pendingWorkItems    count of input tasks
       * @param[in]   inputBuffer         query sequences
       * @param[in]   outputBuffer        container for saving output mappings
       * @param[in]   allReadMappings     container for saving all mappings (optional use depending on filter)
       * @param[in]   totalReadsMapped    counter to track count of reads mapped
       * @param[in]   outstrm             outstream stream object 
       */
      template <typename Vec1, typename Vec2, typename Vec3>
        void processMappingWorkInParallel(int pendingWorkItems, Vec1 &inputBuffer, Vec2 &outputBuffer,
            Vec3 &allReadMappings, seqno_t &totalReadsMapped, std::ofstream &outstrm)
        {
#pragma omp parallel for 
          for(int i=0; i < pendingWorkItems; i++)
          {
            this->mapModule(inputBuffer[i], outputBuffer[i]);

#pragma omp critical
            this->mapModuleHandleOutput(outputBuffer[i], allReadMappings, totalReadsMapped, outstrm);

            //Clear the output container
            outputBuffer[i].reset();
          }
        }

      /**
       * @brief               main mapping function given an input read
       * @details             this function is run in parallel by multiple threads
       * @param[in]   input   input read details
       * @param[out]  output  output object containing the mappings
       */
      void mapModule (MapModuleInput &input, MapModuleOutput &output)
      {
        //save query sequence name
        output.qseqName = input.qseqName;

        if(! param.split)   
        {
          QueryMetaData <MinVec_Type> Q;
          Q.seq = &(input.qseq)[0u];
          Q.len = input.qlen;
          Q.seqCounter = input.seqCounter;

          MappingResultsVector_t l2Mappings;   

          //Map this sequence
          mapSingleQueryFrag(Q, l2Mappings);

          output.readMappings.insert(output.readMappings.end(), l2Mappings.begin(), l2Mappings.end());
        }
        else  //Split read mapping
        {
          int fragmentCount = input.qlen / (param.minMatchLength/2);
          bool mappingReported = false;

          //Map individual short fragments in the read
          for (int i = 0; i < fragmentCount; i++)
          {
            //Prepare fragment sequence object 
            QueryMetaData <MinVec_Type> Q;
            Q.seq = &(input.qseq)[0u] + i * (param.minMatchLength/2);
            Q.len = (param.minMatchLength/2);
            Q.seqCounter = input.seqCounter;

            MappingResultsVector_t l2Mappings;   

            //Map this fragment
            mapSingleQueryFrag(Q, l2Mappings);

            //Adjust query coordinates and length in the reported mapping
            std::for_each(l2Mappings.begin(), l2Mappings.end(), [&](MappingResult &e){ 
                e.queryLen = input.qlen;
                e.queryStartPos = i * (param.minMatchLength/2);
                e.queryEndPos = i * (param.minMatchLength/2) + Q.len - 1;
                });

            output.readMappings.insert(output.readMappings.end(), l2Mappings.begin(), l2Mappings.end());
          }

          //merge
          mergeMappings(output.readMappings);
        }

        //filter mappings best over query sequence axis
        if(param.filterMode == filter::MAP || param.filterMode == filter::ONETOONE)
        {
          skch::Filter::query::filterMappings(output.readMappings);
        }
      }

      /**
       * @brief                       routine to handle mapModule's output of mappings
       * @param[in] output            mapping output object
       * @param[in] allReadMappings   vector to store mappings of all reads (optional use depending on filter)
       * @param[in] totalReadsMapped  counter to track count of reads mapped
       * @param[in] outstrm           outstream stream object 
       */
      template <typename Vec>
        void mapModuleHandleOutput(MapModuleOutput &output, Vec &allReadMappings, seqno_t &totalReadsMapped,
            std::ofstream &outstrm)
        {
          if(output.readMappings.size() > 0)
            totalReadsMapped++;

          if (param.filterMode == filter::ONETOONE)
          {
            //Save for another filtering round
            allReadMappings.insert(allReadMappings.end(), output.readMappings.begin(), output.readMappings.end());
          }
          else
          {  
            //Report mapping
            reportReadMappings(output, outstrm); 
          }

          //Clear mappings in output for re-using object in next iteration
          output.reset();
        }

      /**
       * @brief                   map the parsed query sequence (L1 and L2 mapping)
       * @param[in]   Q           metadata about query sequence
       * @param[in]   outstrm     outstream stream where mappings will be reported
       * @param[out]  l2Mappings  Mapping results in the L2 stage
       */
      template<typename Q_Info, typename VecOut>
        void mapSingleQueryFrag(Q_Info &Q, VecOut &l2Mappings)
        {
#ifdef ENABLE_TIME_PROFILE_L1_L2
          auto t0 = skch::Time::now();
#endif
          //L1 Mapping
          std::vector<L1_candidateLocus_t> l1Mappings; 
          doL1Mapping(Q, l1Mappings);

#ifdef ENABLE_TIME_PROFILE_L1_L2
          std::chrono::duration<double> timeSpentL1 = skch::Time::now() - t0;
          auto t1 = skch::Time::now();
#endif

          //L2 Mapping
          doL2Mapping(Q, l1Mappings, l2Mappings);


#ifdef ENABLE_TIME_PROFILE_L1_L2
          {
            std::chrono::duration<double> timeSpentL2 = skch::Time::now() - t1;
            std::chrono::duration<double> timeSpentMappingFragment = skch::Time::now() - t0;

            std::cerr << Q.seqCounter << " " << Q.len
              << " " << timeSpentL1.count() 
              << " " << timeSpentL2.count()
              << " " << timeSpentMappingFragment.count()
              << "\n";
          }
#endif
        }

      /**
       * @brief       Find candidate regions for a read using level 1 (seed-hits) mapping
       * @details     The count of hits that should occur within a region on the reference is 
       *              determined by the threshold similarity
       *              The resulting start and end target offsets on reference is (are) an 
       *              overestimate of the mapped region. Computing better bounds is left for
       *              the following L2 stage.
       * @param[in]   Q                         query sequence details 
       * @param[out]  l1Mappings                all the read mapping locations
       */
      template <typename Q_Info, typename Vec>
        void doL1Mapping(Q_Info &Q, Vec &l1Mappings)
        {
          //Vector of positions of all the hits 
          std::vector<MinimizerMetaData> seedHitsL1;

          ///1. Compute the minimizers

          CommonFunc::addMinimizers(Q.minimizerTableQuery, Q.seq, Q.len, param.kmerSize, param.windowSize, param.alphabetSize, Q.seqCounter);

#ifdef DEBUG
          std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", minimizer count = " << Q.minimizerTableQuery.size() << "\n";
#endif

          ///2. Find the hits in the reference, pick 's' unique minimizers as seeds, 

          std::sort(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::lessByHash);

          //note : unique preserves the original relative order of elements 
          auto uniqEndIter = std::unique(Q.minimizerTableQuery.begin(), Q.minimizerTableQuery.end(), MinimizerInfo::equalityByHash);

          //This is the sketch size for estimating jaccard
          Q.sketchSize = std::distance(Q.minimizerTableQuery.begin(), uniqEndIter);

          //For invalid query (example : just NNNs), we may be left with 0 sketch size
          //Ignore the query in this case
          if(Q.sketchSize == 0)
            return;

          int totalMinimizersPicked = 0;

          for(auto it = Q.minimizerTableQuery.begin(); it != uniqEndIter; it++)
          {
            //Check if hash value exists in the reference lookup index
            auto seedFind = refSketch.minimizerPosLookupIndex.find(it->hash);

            if(seedFind != refSketch.minimizerPosLookupIndex.end())
            {
              auto hitPositionList = seedFind->second;

              //Save the positions (Ignore high frequency hits)
              if(hitPositionList.size() < refSketch.getFreqThreshold())
              {
                seedHitsL1.insert(seedHitsL1.end(), hitPositionList.begin(), hitPositionList.end());
              }

            }
          }

          int minimumHits = Stat::estimateMinimumHitsRelaxed(Q.sketchSize, param.kmerSize, param.percentageIdentity);

          this->computeL1CandidateRegions(Q, seedHitsL1, minimumHits, l1Mappings);

#ifdef DEBUG
          std::cout << "INFO, skch::Map:doL1Mapping, read id " << Q.seqCounter << ", Count of L1 hits in the reference = " << seedHitsL1.size() << ", minimum hits required for a candidate = " << minimumHits << ", Count of L1 candidate regions = " << l1Mappings.size() << "\n";
#endif

        }

      /**
       * @brief                     Helper function to doL1Mapping()
       * @param[in]   Q             query
       * @param[in]   seedHitsL1    minimizer hits in the reference
       * @param[in]   minimumHits   estimated minimum hits required for significant match
       * @param[out]  l1Mappings    all the read mapping locations
       */
      template <typename Q_Info, typename Vec1, typename Vec2>
        void computeL1CandidateRegions(Q_Info &Q, Vec1 &seedHitsL1, int minimumHits, Vec2 &l1Mappings)
        {
          if(minimumHits < 1)
            minimumHits = 1;

          //Sort all the hit positions
          std::sort(seedHitsL1.begin(), seedHitsL1.end());

          for(auto it = seedHitsL1.begin(); it != seedHitsL1.end(); it++)
          {
            if(std::distance(it, seedHitsL1.end()) >= minimumHits)
            {
              auto it2 = it + minimumHits -1;
              //[it .. it2] are 'minimumHits' consecutive hits 

              //Check if consecutive hits are close enough
              //NOTE: hits may span more than a read length for a valid match, as we keep window positions 
              //      for each minimizer
              if(it2->seqId == it->seqId && it2->wpos - it->wpos < Q.len)
              {
                //Save <1st pos --- 2nd pos>
                L1_candidateLocus_t candidate{it->seqId, 
                    std::max(0, it2->wpos - Q.len + 1), it->wpos};

                //Check if this candidate overlaps with last inserted one
                auto lst = l1Mappings.end(); lst--;

                //match seq_no and see if this candidate begins before last element ends
                if( l1Mappings.size() > 0 
                    && candidate.seqId == lst->seqId 
                    && lst->rangeEndPos >= candidate.rangeStartPos)
                {
                  //Push the end pos of last candidate locus further out
                  lst->rangeEndPos = std::max(candidate.rangeEndPos, lst->rangeEndPos);
                }
                else
                  l1Mappings.push_back(candidate);
              }
            }
          }
        }

      /**
       * @brief                                 Revise L1 candidate regions to more precise locations
       * @param[in]   Q                         query sequence information
       * @param[in]   l1Mappings                candidate regions for query sequence found at L1
       * @param[out]  l2Mappings                Mapping results in the L2 stage
       */
      template <typename Q_Info, typename VecIn, typename VecOut>
        void doL2Mapping(Q_Info &Q, VecIn &l1Mappings, VecOut &l2Mappings)
        {
          ///2. Walk the read over the candidate regions and compute the jaccard similarity with minimum s sketches
          for(auto &candidateLocus: l1Mappings)
          {
            L2_mapLocus_t l2 = {};
            computeL2MappedRegions(Q, candidateLocus, l2);

            //Compute mash distance using calculated jaccard
            float mash_dist = Stat::j2md(1.0 * l2.sharedSketchSize/Q.sketchSize, param.kmerSize);

            //Compute lower bound to mash distance within 90% confidence interval
            float mash_dist_lower_bound = Stat::md_lower_bound(mash_dist, Q.sketchSize, param.kmerSize, 0.75);

            float nucIdentity = 100 * (1 - mash_dist);
            float nucIdentityUpperBound = 100 * (1 - mash_dist_lower_bound);

            //Report the alignment
            if(nucIdentityUpperBound >= param.percentageIdentity)
            {
              MappingResult res;

              //Save the output
              {
                res.queryLen = Q.len;
                res.refStartPos = l2.meanOptimalPos ;
                res.refEndPos = l2.meanOptimalPos + Q.len - 1;
                res.queryStartPos = 0;
                res.queryEndPos = Q.len - 1;
                res.refSeqId = l2.seqId;
                res.querySeqId = Q.seqCounter;
                res.nucIdentity = nucIdentity;
                res.nucIdentityUpperBound = nucIdentityUpperBound;
                res.sketchSize = Q.sketchSize;
                res.conservedSketches = l2.sharedSketchSize;

                //Compute additional statistics -> strand, reference compexity
                {
                  SlideMapper<Q_Info> slidemap(Q);
                  slidemap.insert_ref(l2.optimalStart, l2.optimalEnd);
                  int strandVotes, uniqueRefHashes;
                  slidemap.computeStatistics(strandVotes, uniqueRefHashes);

                  res.strand = strandVotes > 0 ? strnd::FWD : strnd::REV;
                }

                l2Mappings.push_back(res);
              }
            }
          }
        }

      /**
       * @brief                                 Find optimal mapping within an L1 candidate
       * @param[in]   Q                         query sequence information
       * @param[in]   candidateLocus            L1 candidate location
       * @param[out]  l2_out                    L2 mapping inside L1 candidate 
       */
      template <typename Q_Info>
        void computeL2MappedRegions(Q_Info &Q, 
            L1_candidateLocus_t &candidateLocus, 
            L2_mapLocus_t &l2_out)
        {
          //Look up L1 candidate's begin in the index
          MIIter_t firstSuperWindowRangeStart = this->refSketch.searchIndex(candidateLocus.seqId, 
              candidateLocus.rangeStartPos);

          //Count of minimizer windows in a super-window
          offset_t countMinimizerWindows = Q.len - (param.windowSize-1) - (param.kmerSize-1); 

          //Look up the end of the first L2 super-window in the index
          MIIter_t firstSuperWindowRangeEnd = this->refSketch.searchIndex(candidateLocus.seqId, 
              firstSuperWindowRangeStart->wpos + countMinimizerWindows);

          //Look up L1 candidate's end in the index
          MIIter_t lastSuperWindowRangeEnd = this->refSketch.searchIndex(candidateLocus.seqId, 
              candidateLocus.rangeEndPos + Q.len);

          //Define map such that it contains only the query minimizers
          //Used to efficiently compute the jaccard similarity between qry and ref
          SlideMapper<Q_Info> slidemap(Q);

          //Initialize iterator over minimizerIndex
          MIIteratorL2 mi_L2iter( firstSuperWindowRangeStart, firstSuperWindowRangeEnd,
              countMinimizerWindows);

          //Insert all the minimizers in the first 'super-window'
          //  [ mi_L2iter.sw_beg, mi_L2iter.sw_end )
          slidemap.insert_ref(mi_L2iter.sw_beg, mi_L2iter.sw_end);

          auto prev_beg_iter = mi_L2iter.sw_beg;
          auto prev_end_iter = mi_L2iter.sw_end;

          int beginOptimalPos, lastOptimalPos;

          while ( std::distance(mi_L2iter.sw_end, lastSuperWindowRangeEnd) > 0)
          {
            assert( std::distance(mi_L2iter.sw_beg, firstSuperWindowRangeStart) <= 0);
            assert( std::distance(mi_L2iter.sw_end, lastSuperWindowRangeEnd  ) >= 0);

            //Check if the previous first minimizer is out of current range
            if (prev_beg_iter != mi_L2iter.sw_beg)
              slidemap.delete_ref(prev_beg_iter);

            //Check if we have new minimizer in the current range
            if (prev_end_iter != mi_L2iter.sw_end)
              slidemap.insert_ref(prev_end_iter);
          
            //Is this sliding window the best we have so far?
            if (slidemap.sharedSketchElements > l2_out.sharedSketchSize)
            {
              l2_out.sharedSketchSize = slidemap.sharedSketchElements;
              l2_out.optimalStart = mi_L2iter.sw_beg;
              l2_out.optimalEnd = mi_L2iter.sw_end;

              //Save the position
              beginOptimalPos = mi_L2iter.sw_beg->wpos;
              lastOptimalPos = mi_L2iter.sw_beg->wpos;
            }
            else if(slidemap.sharedSketchElements == l2_out.sharedSketchSize)
            {
              //Still save the position
              lastOptimalPos = mi_L2iter.sw_beg->wpos; 
            }

            //Back up the current iterator values
            prev_beg_iter = mi_L2iter.sw_beg;
            prev_end_iter = mi_L2iter.sw_end;

            //Advance the current super-window
            mi_L2iter.next();

          }//End of while loop

          //Save reference sequence id in the mapping output 
          l2_out.seqId = candidateLocus.seqId;
          l2_out.meanOptimalPos = (beginOptimalPos + lastOptimalPos)/2;
        }

      /**
       * @brief                       Merge the consecutive fragment mappings reported in each query 
       * @param[in/out] readMappings  Mappings computed by Mashmap (L2 stage) for a read
       */
      template <typename VecIn>
        void mergeMappings(VecIn &readMappings)
        {
          assert(param.split == true);

          if(readMappings.size() < 2)
            return;

          //length of the read fragment done for split mapping
          auto fragmentLength = param.minMatchLength/2; 

          //Sort the mappings by reference position
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
              return std::tie(a.refSeqId, a.refStartPos, a.queryStartPos) < std::tie(b.refSeqId, b.refStartPos, b.queryStartPos);
              });

          //First assign a unique id to each split mapping in the sorted order
          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            it->splitMappingId = std::distance(readMappings.begin(), it);
          }

          //Start the procedure to identify the chains
          for(auto it = readMappings.begin(); it != readMappings.end(); it++)
          {
            //Which fragment is this wrt. the complete read
            auto currMappingFragno = it->queryStartPos/fragmentLength;

            for(auto it2 = std::next(it); it2 != readMappings.end(); it2++)
            {
              auto thisMappingFragno = it2->queryStartPos / fragmentLength;

              //If this mapping is too far from current mapping being evaluated, stop finding a merge
              if(it2->refSeqId != it->refSeqId || it2->refStartPos - it->refEndPos > 2 * fragmentLength)
                break;

              //If the next mapping is within range, check if it is consecutive query fragment and strand matches
              if( it2->strand == it->strand
                  && thisMappingFragno == currMappingFragno + (it->strand == strnd::FWD ? 1 : -1) )
              {
                it2->splitMappingId = it->splitMappingId;   //merge
                continue;
              }
            }
          }

          //Keep single mapping for each chain and discard others

          //Sort the mappings by post-merge split mapping id
          std::sort(readMappings.begin(), readMappings.end(), [](const MappingResult &a, const MappingResult &b)  
              {
              return a.splitMappingId < b.splitMappingId;
              });

          for(auto it = readMappings.begin(); it != readMappings.end();)
          {
            //Bucket by each chain
            auto it_end = std::find_if(it, readMappings.end(), [&](const MappingResult &e){return e.splitMappingId != it->splitMappingId;} ); 

            //[it -- it_end) represents same chain

            //Incorporate chain information into first mapping

            //compute chain length
            std::for_each(it, it_end, [&](MappingResult &e)
            {
              it->queryStartPos = std::min( it->queryStartPos, e.queryStartPos);
              it->refStartPos = std::min( it->refStartPos, e.refStartPos);

              it->queryEndPos = std::max( it->queryEndPos, e.queryEndPos);
              it->refEndPos = std::max( it->refEndPos, e.refEndPos);
            });

            //Mean identity of all mappings in the chain
            it->nucIdentity = (   std::accumulate(it, it_end, 0.0, 
                                  [](double x, MappingResult &e){ return x + e.nucIdentity; })     )/ std::distance(it, it_end);

            //Discard other mappings of this chain
            std::for_each( std::next(it), it_end, [&](MappingResult &e){ e.discard = 1; });

            //advance the iterator
            it = it_end;
          }

          readMappings.erase( 
              std::remove_if(readMappings.begin(), readMappings.end(), [&](MappingResult &e){ return e.discard == 1; }),
              readMappings.end());
       }

      /**
       * @brief                     Report the final read mappings to output stream
       * @param[in]   output        mapping output i.e. mapping results for a read
       * @param[in]   outstrm       file output stream object
       */
      void reportReadMappings(const MapModuleOutput &output, 
          std::ofstream &outstrm)
      {
        //Print the results
        for(auto &e : output.readMappings)
        {
          assert(e.refSeqId < this->refSketch.metadata.size());

          outstrm  << output.qseqName 
            << " " << e.queryLen 
            << " " << e.queryStartPos
            << " " << e.queryEndPos
            << " " << (e.strand == strnd::FWD ? "+" : "-") 
            << " " << this->refSketch.metadata[e.refSeqId].name
            << " " << this->refSketch.metadata[e.refSeqId].len
            << " " << e.refStartPos 
            << " " << e.refEndPos
            << " " << e.nucIdentity;

          //Print some additional statistics
          outstrm << " " << e.conservedSketches 
            << " " << e.sketchSize << "\n";

          //User defined processing of the results
          if(processMappingResults != nullptr)
            processMappingResults(e);
        }
      }

      /**
       * @brief                         Report the final read mappings to output stream
       * @param[in]   allReadMappings   mapping results for all reads
       * @param[in]   outstrm           file output stream object
       */
      void reportReadMappings(MappingResultsVector_t &allReadMappings, 
          std::ofstream &outstrm)
      {
        //Print the results
        for(auto &e : allReadMappings)
        {
          assert(e.querySeqId < allReadMappings.size());
          assert(e.refSeqId < this->refSketch.metadata.size());

          outstrm  << qmetadata[e.querySeqId].name
            << " " << e.queryLen 
            << " " << e.queryStartPos
            << " " << e.queryEndPos
            << " " << (e.strand == strnd::FWD ? "+" : "-") 
            << " " << this->refSketch.metadata[e.refSeqId].name
            << " " << this->refSketch.metadata[e.refSeqId].len
            << " " << e.refStartPos 
            << " " << e.refEndPos
            << " " << e.nucIdentity;

          //Print some additional statistics
          outstrm << " " << e.conservedSketches 
            << " " << e.sketchSize << "\n";

          //User defined processing of the results
          if(processMappingResults != nullptr)
            processMappingResults(e);
        }
      }

    public:

      /**
       * @brief     An optional utility function to save the 
       *            reported results by the L2 stage into a vector
       */
      static void insertL2ResultsToVec(MappingResultsVector_t &v, const MappingResult &reportedL2Result)
      {
        v.push_back(reportedL2Result);
      }

  };

}

#endif
