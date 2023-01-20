/**
 * @file    commonFunc.hpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMMON_FUNC_HPP
#define COMMON_FUNC_HPP

#include <vector>
#include <map>
#include <algorithm>
#include <deque>
#include <cmath>
#include <fstream>
#include <limits>
#include <numeric>

//Own includes
#include "map/include/map_parameters.hpp"
#include "common/rkmh.hpp"

//External includes
#include "common/murmur3.h"
#include "common/prettyprint.hpp"

#include "assert.hpp"

namespace skch {
    /**
     * @namespace skch::CommonFunc
     * @brief     Implements frequently used common functions
     */
    namespace CommonFunc {
        //seed for murmerhash
        const int seed = 42;

        // Pivot to keep track of sketch border
        template <typename I>
        struct Pivot {
            I p;
            int64_t rank;
        };

        /**
         * @brief   reverse complement of kmer (borrowed from mash)
         * @note    assumes dest is pre-allocated
         */
        inline void reverseComplement(const char *src, char *dest, int length) {
            for (int i = 0; i < length; i++) {
                char base = src[i];

                switch (base) {
                    case 'A':
                        base = 'T';
                        break;
                    case 'C':
                        base = 'G';
                        break;
                    case 'G':
                        base = 'C';
                        break;
                    case 'T':
                        base = 'A';
                        break;
                    default:
                        break;
                }

                dest[length - i - 1] = base;
            }
        }

        /**
     * @brief               convert DNA or AA alphabets to upper case, converting non-canonical DNA bases to N
     * @param[in]   seq     pointer to input sequence
     * @param[in]   len     length of input sequence
     */
        inline void makeUpperCaseAndValidDNA(char *seq, offset_t len) {
            for (int i = 0; i < len; i++) {
                if (seq[i] > 96 && seq[i] < 123) {
                    seq[i] -= 32;
                }

                if (rkmh::valid_dna[seq[i]]) {
                    seq[i] = 'N';
                }
            }
        }

//        /**
//       * @brief               convert DNA or AA alphabets to upper case
//       * @param[in]   seq     pointer to input sequence
//       * @param[in]   len     length of input sequence
//       */
//        inline void makeUpperCase(char *seq, offset_t len) {
//            for (int i = 0; i < len; i++) {
//                if (seq[i] > 96 && seq[i] < 123) {
//                    seq[i] -= 32;
//                }
//            }
//        }
//
//        /**
//         * @brief               convert non-canonical DNA bases to N
//         * @param[in]   seq     pointer to input sequence
//         * @param[in]   len     length of input sequence
//         */
//        inline void makeValidDNA(char *seq, offset_t len) {
//            for (int i = 0; i < len; i++) {
//                if (rkmh::valid_dna[seq[i]]) {
//                    seq[i] = 'N';
//                }
//            }
//        }

        /**
         * @brief   hashing kmer string (borrowed from mash)
         */
        inline hash_t getHash(const char *seq, int length) {
            char data[16];
            MurmurHash3_x64_128(seq, length, seed, data);

            hash_t hash;

            hash = *((hash_t *) data);

            return hash;
        }

        /**
         * @brief		takes hash value of kmer and adjusts it based on kmer's weight
         *					this value will determine its order for minimizer selection
         * @details	this is inspired from Chum et al.'s min-Hash and tf-idf weighting
         */
//        static inline double applyWeight(char* kmer, int kmer_size, hash_t kmer_hash, const std::unordered_set<std::string>& high_freq_kmers) {
//            double x = kmer_hash * 1.0 / UINT32_MAX;  //bring it within [0, 1]
//            //assert (x >= 0.0 && x <= 1.0);
//
//            std::string kmer_str(kmer, kmer_size);
//            if (high_freq_kmers.count(kmer_str) > 0) {
//                /* downweigting by a factor of 8 */
//                /* further aggressive downweigting may affect accuracy */
//                double p2 = x*x;
//                double p4 = p2 * p2;
//                return - 1.0 * (p4 * p4);
//            }
//            return -1.0 * x;
//
//            //range of returned value is between [-1,0]
//            //we avoid adding one for better double precision
//        }


        /**
         * @brief       Compute the minimum s kmers for a string.
         * @param[out]  minmerIndex     container storing sketched Kmers 
         * @param[in]   seq                 pointer to input sequence
         * @param[in]   len                 length of input sequence
         * @param[in]   kmerSize
         * @param[in]   s                   sketch size. 
         * @param[in]   seqCounter          current sequence number, used while saving the position of minimizer
         */
        template <typename T>
          inline void sketchSequence(
              std::vector<T> &minmerIndex, 
              char* seq, 
              offset_t len,
              int kmerSize, 
              int alphabetSize,
              int sketchSize,
              seqno_t seqCounter)
        {
          makeUpperCaseAndValidDNA(seq, len);

          //Compute reverse complement of seq
          std::unique_ptr<char[]> seqRev(new char[len]);
          //char* seqRev = new char[len];

          if(alphabetSize == 4) //not protein
            CommonFunc::reverseComplement(seq, seqRev.get(), len);

          // TODO cleanup
          std::unordered_map<hash_t, MinmerInfo> sketched_vals;
          std::vector<hash_t> sketched_heap;

          for(offset_t i = 0; i < len - kmerSize + 1; i++)
          {
            //Hash kmers
            hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize); 
            hash_t hashBwd;

            if(alphabetSize == 4)
              hashBwd = CommonFunc::getHash(seqRev.get() + len - i - kmerSize, kmerSize);
            else  //proteins
              hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

            //Consider non-symmetric kmers only
            if(hashBwd != hashFwd)
            {
              //Take minimum value of kmer and its reverse complement
              hash_t currentKmer = std::min(hashFwd, hashBwd);

              //Check the strand of this minimizer hash value
              auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

              ////std::cout << seqCounter << "\t" << std::string(seq + i, kmerSize) << "-->" <<  currentKmer << std::endl;
              if (sketched_vals.empty() || sketched_vals.find(currentKmer) == sketched_vals.end()) 
              {

                // Add current hash to heap
                if (sketched_vals.size() < sketchSize || currentKmer < sketched_heap[0])  
                {
                    sketched_vals[currentKmer] = MinmerInfo{currentKmer, seqCounter, i, i, currentStrand};
                    sketched_heap.push_back(currentKmer);
                    std::push_heap(sketched_heap.begin(), sketched_heap.end());
                }

                // Remove one if too large
                if (sketched_vals.size() > sketchSize) 
                {
                    sketched_vals.erase(sketched_heap[0]);
                    std::pop_heap(sketched_heap.begin(), sketched_heap.end());
                    sketched_heap.pop_back();
                }
              } 
              else 
              {
                // TODO these sketched values might never be useful, might save memory by deleting
                // extend the length of the window
                sketched_vals[currentKmer].wpos_end = i;
                sketched_vals[currentKmer].strand += currentStrand == strnd::FWD ? 1 : -1;
              }
            }
          }
          //DEBUG_ASSERT(sketched_vals.size() <= sketchSize);
          minmerIndex.reserve(sketched_vals.size());
          std::for_each(sketched_vals.begin(), sketched_vals.end(),
              [&minmerIndex](auto& pair) {
              pair.second.strand = pair.second.strand > 0 ? strnd::FWD : (pair.second.strand == 0 ? strnd::AMBIG : strnd::REV);
              minmerIndex.push_back(pair.second);
          });

          return;
        }
        

        /**
         * @brief       Compute winnowed minmers from a given sequence and add to the index
         * @param[out]  minmerIndex  table storing minmers and their position as we compute them
         * @param[in]   seq             pointer to input sequence
         * @param[in]   len             length of input sequence
         * @param[in]   kmerSize
         * @param[in]   windowSize
         * @param[in]   sketchSize      sketch size. 
         * @param[in]   seqCounter      current sequence number, used while saving the position of minimizer
         */
        template <typename T>
          inline void addMinmers(std::vector<T> &minmerIndex, 
              char* seq, offset_t len,
              int kmerSize, 
              int windowSize,
              int alphabetSize,
              int sketchSize,
              seqno_t seqCounter)
          {
            /**
             * Double-ended queue (saves minimum at front end)
             * Saves pair of the minimizer and the position of hashed kmer in the sequence
             * Position of kmer is required to discard kmers that fall out of current window
             */
            std::deque< std::tuple<hash_t, strand_t, offset_t> > Q;
            using windowMap_t = std::map<hash_t, std::pair<MinmerInfo, uint64_t>>;
            windowMap_t sortedWindow;
            Pivot<typename windowMap_t::iterator>  piv = {sortedWindow.begin(), 0};

            makeUpperCaseAndValidDNA(seq, len);

            //Compute reverse complement of seq
            std::unique_ptr<char[]> seqRev(new char[kmerSize]);

            //if(alphabetSize == 4) //not protein
              //CommonFunc::reverseComplement(seq, seqRev.get(), len);

            for(offset_t i = 0; i < len - kmerSize + 1; i++)
            {
              //The serial number of current sliding window
              //First valid window appears when i = windowSize - 1
              offset_t currentWindowId = i + kmerSize - windowSize;

              if (currentWindowId == 0) 
              {
                uint64_t rank = 1;
                auto iter = sortedWindow.begin();
                while (iter != sortedWindow.end() && rank <= sketchSize) 
                {
                  iter->second.first.wpos = currentWindowId; 
                  std::advance(iter, 1);
                  rank += 1;
                }
              }

              //Hash kmers
              hash_t hashFwd = CommonFunc::getHash(seq + i, kmerSize); 
              hash_t hashBwd;

              if(alphabetSize == 4) 
              {
                CommonFunc::reverseComplement(seq + i, seqRev.get(), kmerSize);
                hashBwd = CommonFunc::getHash(seqRev.get(), kmerSize);
              }
              else  //proteins
                hashBwd = std::numeric_limits<hash_t>::max();   //Pick a dummy high value so that it is ignored later

              //Take minimum value of kmer and its reverse complement
              hash_t currentKmer = std::min(hashFwd, hashBwd);

              //Check the strand of this minimizer hash value
              auto currentStrand = hashFwd < hashBwd ? strnd::FWD : strnd::REV;

              //If front minimum is not in the current window, remove it
              if (!Q.empty() && std::get<2>(Q.front()) <  currentWindowId) 
              {
                const auto [leaving_hash, leaving_strand, _] = Q.front();

                // Check if we've deleted the hash already
                if (sortedWindow.find(leaving_hash) != sortedWindow.end()) 
                {
                  // If the hash that is getting popped off is still in the window and it is now leaving the window 
                  // wpos != -1 and wpos_end == -1 --> still in window
                  if (sortedWindow[leaving_hash].first.wpos != -1 and sortedWindow[leaving_hash].first.wpos_end == -1 && sortedWindow[leaving_hash].second == 1) 
                  {
                    sortedWindow[leaving_hash].first.wpos_end = currentWindowId;
                    minmerIndex.push_back(sortedWindow[leaving_hash].first);
                  } 

                  // Remove hash
                  sortedWindow[leaving_hash].second -= 1;
                  if (sortedWindow[leaving_hash].second == 0) 
                  {
                    if (leaving_hash == piv.p->first) 
                    {
                      std::advance(piv.p, 1);
                    }
                    else if (leaving_hash < piv.p->first) 
                    {
                      // Kicking out a sketched element
                      if (sortedWindow.size() >= sketchSize + 1) 
                      {
                          std::advance(piv.p, 1);
                      }
                    }
                    sortedWindow.erase(leaving_hash);
                  } 
                  else 
                  {
                    // Not removing hash, but need to adjust the strand
                    if ((sortedWindow[leaving_hash].first.strand == 0 || sortedWindow[leaving_hash].first.strand - leaving_strand == 0)
                        && leaving_hash < piv.p->first) 
                    {
                      sortedWindow[leaving_hash].first.wpos_end = currentWindowId;
                      minmerIndex.push_back(sortedWindow[leaving_hash].first);
                      sortedWindow[leaving_hash].first.wpos = currentWindowId;
                      sortedWindow[leaving_hash].first.wpos_end = -1;
                    }
                    sortedWindow[leaving_hash].first.strand -= leaving_strand;
                  }
                }
                Q.pop_front();
              }

              //Consider non-symmetric kmers only
              if(hashBwd != hashFwd)
              {
                // Add current hash to window
                Q.push_back(std::make_tuple(currentKmer, currentStrand, i)); 
                if (sortedWindow[currentKmer].second == 0) {
                  auto mi = MinmerInfo{currentKmer, seqCounter, -1, -1, currentStrand};
                  sortedWindow[currentKmer].first = mi;

                  if (sortedWindow.size() >= sketchSize + 2 && currentKmer < piv.p->first) {
                      piv.p--;
                  }
                  if (sortedWindow.size() <= sketchSize) {
                      piv.p = sortedWindow.end();
                  } else if (sortedWindow.size() == sketchSize + 1) {
                      piv.p = std::prev(sortedWindow.end());
                  }
                } else {
                  if ((sortedWindow[currentKmer].first.strand + currentStrand == 0 || sortedWindow[currentKmer].first.strand == 0) 
                      && currentKmer < piv.p->first) {
                    sortedWindow[currentKmer].first.wpos_end = currentWindowId;
                    minmerIndex.push_back(sortedWindow[currentKmer].first);
                    sortedWindow[currentKmer].first.wpos = currentWindowId;
                    sortedWindow[currentKmer].first.wpos_end = -1;
                  }
                  sortedWindow[currentKmer].first.strand += currentStrand;
                }
                sortedWindow[currentKmer].second += 1;
              }

              //Select the minimizer from Q and put into index
              if(currentWindowId >= 0)
              {

                // Does the new kmer belong in the sketch?
                if (hashBwd != hashFwd                                                  // Non-symmetric 
                    && ((piv.p == sortedWindow.end()) || (currentKmer < piv.p->first))  // Belongs in sketch
                    && sortedWindow[currentKmer].first.wpos == -1)  // Haven't seen it in the window yet
                {                                     
                  sortedWindow[currentKmer].first.wpos = currentWindowId;
                }

                // Did we incorporate a previously hashed kmer into the sketch?
                if (sortedWindow.size() != 0) 
                {
                  auto& sth_mi = std::prev(piv.p)->second.first;
                  if (sth_mi.wpos == -1) 
                  {
                    sth_mi.wpos = currentWindowId;
                  }
                }

                // Did we kick a minmer into non-sketch territory?
                if (piv.p != sortedWindow.end()) 
                {
                  auto& splus1th_mi = piv.p->second.first;
                  if (splus1th_mi.wpos != -1) 
                  {
                    splus1th_mi.wpos_end = currentWindowId;
                    minmerIndex.push_back(MinmerInfo(splus1th_mi));
                    splus1th_mi.wpos = -1;
                    splus1th_mi.wpos_end = -1;
                  }
                }
#ifdef DEBUG
                // Brute force ensure minmer validity
                //DEBUG_ASSERT(std::distance(sortedWindow.begin(), piv.p) == std::min<int>(sortedWindow.size(), sketchSize), seqCounter, currentWindowId, i);
                //DEBUG_ASSERT(piv.p == sortedWindow.end() || (piv.p->second.first.wpos == -1 && piv.p->second.first.wpos_end == -1));
                //DEBUG_ASSERT((sortedWindow.size() == 0 || currentWindowId < 0) || (std::prev(piv.p)->second.first.wpos != -1 && std::prev(piv.p)->second.first.wpos_end == -1));
                //for (auto it = sortedWindow.begin(); it != sortedWindow.end(); it++) {
                  //if (piv.p == sortedWindow.end() || it->first < piv.p->first) {
                    //DEBUG_ASSERT(it->second.first.wpos != -1, it->second.first);
                    //DEBUG_ASSERT(it->second.first.wpos_end == -1);
                  //} else {
                    //DEBUG_ASSERT(it->second.first.wpos == -1, it->second.first, currentWindowId);
                    //DEBUG_ASSERT(it->second.first.wpos_end == -1);
                  //}
                //}
#endif
              }
              else 
              {
                if (hashBwd != hashFwd && sortedWindow[currentKmer].second == 1) 
                {
                  // Seeing kmer for the first time
                  if (sortedWindow.size() < sketchSize + 1) 
                  {
                    piv.p = sortedWindow.end();
                  } 
                  else if (sortedWindow.size() == sketchSize + 1) 
                  {
                    piv.p = std::prev(sortedWindow.end());
                  }
                }
              }

              //DEBUG_ASSERT(sortedWindow.size() == 0 
                      //|| std::distance(sortedWindow.begin(), piv.p) == std::min<int>(sortedWindow.size(), sketchSize), 
                      //seqCounter, currentWindowId, i);
              //DEBUG_ASSERT(piv.p == sortedWindow.end() || (piv.p->second.first.wpos == -1 && piv.p->second.first.wpos_end == -1));
              //DEBUG_ASSERT(((sortedWindow.size() == 0 || currentWindowId < 0) || (std::prev(piv.p)->second.first.wpos != -1 && std::prev(piv.p)->second.first.wpos_end == -1)));
            }

            uint64_t rank = 1;
            auto iter = sortedWindow.begin();
            while (iter != sortedWindow.end() && rank <= sketchSize) 
            {
              if (iter->second.first.wpos != -1) 
              {
                iter->second.first.wpos_end = len - kmerSize + 1;
                minmerIndex.push_back(iter->second.first);
              }
              std::advance(iter, 1);
              rank += 1;
            }

            // TODO Not sure why these are occuring but they are a bug
            minmerIndex.erase(
                std::remove_if(
                  minmerIndex.begin(), 
                  minmerIndex.end(), 
                  [](auto& mi) { return mi.wpos < 0 || mi.wpos_end < 0; }),
                minmerIndex.end());


            // Split up windows longer than windowSize into chunks of windowSize or less
            std::vector<MinmerInfo> chunkedMIs;
            std::for_each(minmerIndex.begin(), minmerIndex.end(), [&chunkedMIs, windowSize] (auto& mi) {
              mi.strand = mi.strand < 0 ? (mi.strand == 0 ? strnd::AMBIG : strnd::REV) : strnd::FWD;
              if (mi.wpos_end > mi.wpos + windowSize) {
                for (int chunk = 0; chunk < std::ceil(float(mi.wpos_end - mi.wpos) / float(windowSize)); chunk++) {
                  chunkedMIs.push_back(
                    MinmerInfo{
                      mi.hash, 
                      mi.seqId, 
                      mi.wpos + chunk*windowSize, 
                      std::min(mi.wpos + chunk*windowSize + windowSize, mi.wpos_end),
                      mi.strand
                    } 
                  );
                }
              }
            });
            minmerIndex.erase(
                std::remove_if(
                  minmerIndex.begin(), 
                  minmerIndex.end(), 
                  [windowSize](auto& mi) { return mi.wpos_end - mi.wpos > windowSize; }),
                minmerIndex.end());
            minmerIndex.insert(minmerIndex.end(), chunkedMIs.begin(), chunkedMIs.end());

            // Sort the index based on start position
            std::sort(minmerIndex.begin(), minmerIndex.end(), [](auto& l, auto& r) {return l.wpos < r.wpos;});

            // No duplicate windows
            // TODO These should not be occurring. They happen rarely, so just deleting them for now
            // but need to fix eventually 
            minmerIndex.erase(
                std::unique(
                  minmerIndex.begin(), 
                  minmerIndex.end(), 
                  [](auto& l, auto& r) { return (l.wpos == r.wpos) && (l.hash == r.hash); }),
                minmerIndex.end());

#ifdef DEBUG
            ////std::cout << "INFO, skch::CommonFunc::addMinimizers, inserted minimizers for sequence id = " << seqCounter << "\n";
            ////std::cout << "INFO, skch::CommonFunc::addMinimizers, length of sequence  = " << len << "\n";
            //DEBUG_ASSERT(std::all_of(minmerIndex.begin(), minmerIndex.end(), [](auto& mi) {return mi.wpos >= 0;}));
            //DEBUG_ASSERT(std::all_of(minmerIndex.begin(), minmerIndex.end(), [](auto& mi) {return mi.wpos_end >= 0;}));
            //std::vector<MinmerInfo> endpos_heap;
            //auto heap_cmp = [](auto& l, auto& r) {return l.wpos_end >= r.wpos_end;};
            //for (auto& mi : minmerIndex) {
              //while (!endpos_heap.empty() && endpos_heap.front().wpos_end <= mi.wpos) {
                //std::pop_heap(endpos_heap.begin(), endpos_heap.end(), heap_cmp); 
                //endpos_heap.pop_back();
              //}
              //endpos_heap.push_back(mi);
              //std::push_heap(endpos_heap.begin(), endpos_heap.end(), heap_cmp);
              //DEBUG_ASSERT(endpos_heap.size() <= sketchSize);
            //}
#endif
          }

        /**
          * @brief           Functor for comparing tuples by single index layer
          * @tparam layer    Tuple's index which is used for comparison
          * @tparam op       comparator, default as std::less
          */
        template<size_t layer, template<typename> class op = std::less>
        struct TpleComp {
            //Compare two tuples using their values
            template<typename T>
            bool operator()(T const &t1, T const &t2) {
                return op<typename std::tuple_element<layer, T>::type>()(std::get<layer>(t1), std::get<layer>(t2));
            }
        };

        /**
         * @brief                   computes the total size of reference in bytes
         * @param[in] refSequences  vector of reference files
         * @return                  total size
         */
        inline uint64_t getReferenceSize(const std::vector<std::string> &refSequences) {
            uint64_t count = 0;

            for (auto &f : refSequences) {
                //Open the file as binary, and set the position to end
                std::ifstream in(f, std::ifstream::ate | std::ifstream::binary);

                //the position of the current character
                count += (uint64_t) (in.tellg());
            }

            return count;
        }

        // Splitting
        template<typename Out>
        void split(const std::string &s, char delim, Out result) {
            std::stringstream ss(s);
            std::string item;
            while (std::getline(ss, item, delim)) {
                *(result++) = item;
            }
        }

        std::vector<std::string> split(const std::string &s, char delim) {
            std::vector<std::string> elems;
            split(s, delim, std::back_inserter(elems));
            return elems;
        }
    }
}

#endif
