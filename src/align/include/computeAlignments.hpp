/**
 * @file    computeAlignments.hpp
 * @brief   logic for generating alignments when given mashmap 
 *          mappings as input
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef COMPUTE_ALIGNMENTS_HPP 
#define COMPUTE_ALIGNMENTS_HPP

#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <zlib.h>
#include <cassert>
#include <thread>

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/kseq.h"
#include "common/edlib.h"
#include "common/atomic_queue/atomic_queue.h"

KSEQ_INIT(gzFile, gzread)

namespace align
{

  long double float2phred(long double prob) {
    if (prob == 1)
        return 255;  // guards against "-0"
    long double p = -10 * (long double) log10(prob);
    if (p < 0 || p > 255) // int overflow guard
        return 255;
    else
        return p;
  }

  struct seq_record_t {
      MappingBoundaryRow currentRecord;
      std::string mappingRecordLine;
      std::string qSequence;
      seq_record_t(const MappingBoundaryRow& c, const std::string& r, const std::string& q)
          : currentRecord(c)
          , mappingRecordLine(r)
          , qSequence(q)
          { }
  };
  // load into this
  typedef atomic_queue::AtomicQueue<seq_record_t*, 2 << 16> seq_atomic_queue_t;
  // results into this, write out
  typedef atomic_queue::AtomicQueue<std::string*, 2 << 16> paf_atomic_queue_t;

  /**
   * @class     align::Aligner
   * @brief     compute alignments and generate sam output
   *            from mashmap mappings
   */
  class Aligner
  {
    private:

      //algorithm parameters
      const align::Parameters &param;

      refSequenceMap_t refSequences;

    public:

      /**
       * @brief                 constructor, also reads reference sequences
       * @param[in] p           algorithm parameters
       */
      Aligner(const align::Parameters &p) :
        param(p)
      {
        this->getRefSequences();
      }

      /**
       * @brief                 compute aligments
       */
      void compute()
      {
        this->computeAlignments();
      }

    private:

      /**
       * @brief                 parse and save all the reference sequences
       */
      void getRefSequences()
      {
        for(const auto &fileName : param.refSequences)
        {

#ifdef DEBUG
          std::cout << "INFO, align::Aligner::getRefSequences, parsing reference sequences in file " << fileName << std::endl;
#endif

          //Open the file using kseq
          FILE *file = fopen(fileName.c_str(), "r");
          gzFile fp = gzdopen(fileno(file), "r");
          kseq_t *seq = kseq_init(fp);


          //size of sequence
          skch::offset_t len;

          while ((len = kseq_read(seq)) >= 0) 
          {
             std::string seqId = seq->name.s;

             skch::CommonFunc::makeUpperCase(seq->seq.s, len);
             std::string sequence = seq->seq.s;

             //seqId shouldn't already exist in our table
             assert(this->refSequences.count(seqId) == 0);

             refSequences.emplace(seqId, sequence);
          }

          kseq_destroy(seq);  
          gzclose(fp); //close the file handler 
          fclose(file);
        }
      }

      char* alignmentToCigar(const unsigned char* const alignment,
                             const int alignmentLength,
                             uint64_t& refAlignedLength,
                             uint64_t& qAlignedLength,
                             uint64_t& matches,
                             uint64_t& mismatches,
                             uint64_t& insertions,
                             uint64_t& deletions,
                             uint64_t& softclips) {

          // Maps move code from alignment to char in cigar.
          //                        0    1    2    3
          char moveCodeToChar[] = {'=', 'I', 'D', 'X'};

          vector<char>* cigar = new vector<char>();
          char lastMove = 0;  // Char of last move. 0 if there was no previous move.
          int numOfSameMoves = 0;
          for (int i = 0; i <= alignmentLength; i++) {
              // if new sequence of same moves started
              if (i == alignmentLength || (moveCodeToChar[alignment[i]] != lastMove && lastMove != 0)) {
                  // calculate matches, mismatches, insertions, deletions
                  switch (lastMove) {
                  case 'I':
                      // assume that starting and ending insertions are softclips
                      if (i == alignmentLength || cigar->empty()) {
                          softclips += numOfSameMoves;
                      } else {
                          insertions += numOfSameMoves;
                      }
                      qAlignedLength += numOfSameMoves;
                      break;
                  case '=':
                      matches += numOfSameMoves;
                      qAlignedLength += numOfSameMoves;
                      refAlignedLength += numOfSameMoves;
                      break;
                  case 'X':
                      mismatches += numOfSameMoves;
                      qAlignedLength += numOfSameMoves;
                      refAlignedLength += numOfSameMoves;
                      break;
                  case 'D':
                      deletions += numOfSameMoves;
                      refAlignedLength += numOfSameMoves;
                      break;
                  default:
                      break;
                  }
                  
                  // Write number of moves to cigar string.
                  int numDigits = 0;
                  for (; numOfSameMoves; numOfSameMoves /= 10) {
                      cigar->push_back('0' + numOfSameMoves % 10);
                      numDigits++;
                  }
                  reverse(cigar->end() - numDigits, cigar->end());
                  // Write code of move to cigar string.
                  cigar->push_back(lastMove);
                  // If not at the end, start new sequence of moves.
                  if (i < alignmentLength) {
                      // Check if alignment has valid values.
                      if (alignment[i] > 3) {
                          delete cigar;
                          return 0;
                      }
                      numOfSameMoves = 0;
                  }
              }
              if (i < alignmentLength) {
                  lastMove = moveCodeToChar[alignment[i]];
                  numOfSameMoves++;
              }
          }
          cigar->push_back(0);  // Null character termination.
          char* cigar_ = (char*) malloc(cigar->size() * sizeof(char));
          memcpy(cigar_, &(*cigar)[0], cigar->size() * sizeof(char));
          delete cigar;

          return cigar_;
      }
      
      /**
       * @brief                 parse query sequences and mashmap mappings
       *                        to compute sequence alignments
       */
      void computeAlignments()
      {

          // input atomic queue
          seq_atomic_queue_t seq_queue;
          // output atomic queue
          paf_atomic_queue_t paf_queue;
          // flag when we're done reading
          std::atomic<bool> reader_done;
          reader_done.store(false);

          auto& nthreads = param.threads;

          // atomics to record if we're working or not
          std::vector<std::atomic<bool>> working(nthreads);
          for (auto& w : working) {
              w.store(true);
          }

          // reader picks up candidate alignments from input
          auto reader_thread =
              [&](void) {
                  //Parse query sequences
                  for(const auto &fileName : param.querySequences)
                  {
#ifdef DEBUG
                      std::cout << "INFO, align::Aligner::computeAlignments, parsing query sequences in file " << fileName << std::endl;
#endif

                      //Open the file using kseq
                      FILE *file = fopen(fileName.c_str(), "r");
                      gzFile fp = gzdopen(fileno(file), "r");
                      kseq_t *seq = kseq_init(fp);

                      //Open mashmap output file
                      std::ifstream mappingListStream(param.mashmapPafFile);
                      std::string mappingRecordLine;
                      MappingBoundaryRow currentRecord;

                      //size of sequence
                      skch::offset_t len;

                      while ((len = kseq_read(seq)) >= 0) 
                      {
                          std::string qSeqId = seq->name.s;

                          skch::CommonFunc::makeUpperCase(seq->seq.s, len);
                          // todo maybe this should change to some kind of unique pointer?
                          // something where we can GC it when we're done aligning to it
                          std::string qSequence = seq->seq.s;

                          //Check if all mapping records are processed already
                          if( mappingListStream.eof() )
                              break;

                          //Read first record from mashmap output file during first iteration
                          if(mappingRecordLine.empty())       
                          {
                              std::getline(mappingListStream, mappingRecordLine);
                          }

                          this->parseMashmapRow(mappingRecordLine, currentRecord);

                          //Check if mapping query id matches current query sequence id
                          if(currentRecord.qId != qSeqId)
                          {
                              //Continue to read the next query sequence
                              continue;
                          }
                          else
                          {
                              auto q = new seq_record_t(currentRecord, mappingRecordLine, qSequence);
                              seq_queue.push(q);
                          }

                          //Check if more mappings have same query sequence id
                          while(std::getline(mappingListStream, mappingRecordLine))
                          {
                              this->parseMashmapRow(mappingRecordLine, currentRecord);

                              if(currentRecord.qId != qSeqId)
                              {
                                  //Break the inner loop to read query sequence
                                  break;
                              }
                              else
                              {
                                  auto q = new seq_record_t(currentRecord, mappingRecordLine, qSequence);
                                  seq_queue.push(q);
                              }
                          }
                      }

                      mappingListStream.close();

                      kseq_destroy(seq);  
                      gzclose(fp); //close the file handler 
                      fclose(file);
                  }
                  reader_done.store(true);
              };

          // helper to check if we're still aligning
          auto still_working =
              [&](const std::vector<std::atomic<bool>>& working) {
                  bool ongoing = false;
                  for (auto& w : working) {
                      ongoing = ongoing || w.load();
                  }
                  return ongoing;
              };

          // writer, picks output from queue and writes it to our output stream
          std::ofstream outstrm(param.samOutputFile);
          auto writer_thread =
              [&](void) {
                  while (true) {
                      std::string* paf_lines = nullptr;
                      if (!paf_queue.try_pop(paf_lines)
                          && !still_working(working)) {
                          break;
                      } else if (paf_lines != nullptr) {
                          outstrm << *paf_lines;
                          delete paf_lines;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
              };

          // worker, takes candidate alignments and runs edlib alignment on them
          auto worker_thread = 
              [&](uint64_t tid,
                  std::atomic<bool>& is_working) {
                  is_working.store(true);
                  while (true) {
                      seq_record_t* rec = nullptr;
                      if (!seq_queue.try_pop(rec)
                          && reader_done.load()) {
                          break;
                      } else if (rec != nullptr) {
                          std::string* paf_rec
                              = new std::string(
                                  doAlignment(rec->currentRecord,
                                              rec->mappingRecordLine,
                                              rec->qSequence));
                          if (paf_rec->size()) {
                              paf_queue.push(paf_rec);
                          } else {
                              delete paf_rec;
                          }
                          delete rec;
                      } else {
                          std::this_thread::sleep_for(100ns);
                      }
                  }
                  is_working.store(false);
              };

          // launch reader
          std::thread reader(reader_thread);
          // launch writer
          std::thread writer(writer_thread);
          // launch workers
          std::vector<std::thread> workers; workers.reserve(nthreads);
          for (uint64_t t = 0; t < nthreads; ++t) {
              workers.emplace_back(worker_thread,
                                   t,
                                   std::ref(working[t]));
          }

          // wait for reader and workers to complete
          reader.join();
          for (auto& worker : workers) {
              worker.join();
          }
          // and finally the writer
          writer.join();

      }

      /**
       * @brief                         parse mashmap row sequence 
       * @param[in]   mappingRecordLine
       * @param[out]  currentRecord
       */
      inline void parseMashmapRow(const std::string &mappingRecordLine, MappingBoundaryRow &currentRecord)
      {
        std::stringstream ss(mappingRecordLine); // Insert the string into a stream
        std::string word; // Have a buffer string

        vector<std::string> tokens; // Create vector to hold our words

        while (ss >> word)
          tokens.push_back(word);

        //We expect and need at least these many values in a mashmap mapping
        assert(tokens.size() >= 9);

        //Save words into currentRecord
        {
          currentRecord.qId = tokens[0];
          currentRecord.qStartPos = std::stoi(tokens[2]);
          currentRecord.qEndPos = std::stoi(tokens[3]);
          currentRecord.strand = (tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV);
          currentRecord.refId = tokens[5];
          currentRecord.rStartPos = std::stoi(tokens[7]);
          currentRecord.rEndPos = std::stoi(tokens[8]);
        }
      }

      /**
       * @brief                           compute alignment using edlib 
       * @param[in]   currentRecord       mashmap mapping parsed information
       * @param[in]   mappingRecordLine   mashmap mapping output raw string
       * @param[in]   qSequence           query sequence
       * @param[in]   outstrm             output stream
       */
      std::string doAlignment(MappingBoundaryRow &currentRecord, const std::string &mappingRecordLine, const std::string &qSequence)
      {

#ifdef DEBUG
        std::cout << "INFO, align::Aligner::doAlignment, aligning mashmap record: " << mappingRecordLine << std::endl;
#endif

        //Define reference substring for this mapping
        const std::string &refId = currentRecord.refId;
        const char* refRegion = this->refSequences[refId].c_str();
        const auto& refSize = this->refSequences[refId].size();
        //currentRecord.rStartPos = std::max((int64_t)0, (int64_t)currentRecord.rStartPos - 1000);
        //currentRecord.rEndPos = std::min((int64_t)refSize, (int64_t)currentRecord.rEndPos + 1000);
        refRegion += currentRecord.rStartPos;
        skch::offset_t refLen = currentRecord.rEndPos - currentRecord.rStartPos + 1;
        assert(refLen <= refSize);

        //Define query substring for this mapping
        const char* queryRegion = qSequence.c_str();  //initially point to beginning
        const auto& querySize = qSequence.size();
        //currentRecord.qStartPos = std::max((int64_t)0, (int64_t)currentRecord.qStartPos - 1000);
        //currentRecord.qEndPos = std::min((int64_t)querySize, (int64_t)currentRecord.qEndPos + 1000);
        skch::offset_t queryLen = currentRecord.qEndPos - currentRecord.qStartPos + 1;
        queryRegion += currentRecord.qStartPos;

        char* queryRegionStrand = new char[queryLen];

        if(currentRecord.strand == skch::strnd::FWD) {
          strncpy(queryRegionStrand, queryRegion, queryLen);    //Copy the same string
        } else {
          skch::CommonFunc::reverseComplement(queryRegion, queryRegionStrand, queryLen); //Reverse complement
        }

        assert(queryLen <= querySize);

        //Compute alignment
        auto t0 = skch::Time::now();

        //for defining size of edlib's band during alignment
        int editDistanceLimit;

        if(param.percentageIdentity == 0)
          editDistanceLimit = -1;   //not bounded
        else
          editDistanceLimit = (int)((1 - param.percentageIdentity/100) * queryLen);

        if(param.bandwidth > 0)
          editDistanceLimit = std::min(editDistanceLimit, param.bandwidth);

#ifdef DEBUG
        std::cout << "INFO, align::Aligner::doAlignment, edlib execution starting, query region length = " << queryLen
          << ", reference region length= " << refLen << ", edit distance limit= " << editDistanceLimit << std::endl; 
#endif

        EdlibAlignResult result = edlibAlign(
            queryRegionStrand, queryLen, refRegion, refLen,
            edlibNewAlignConfig(editDistanceLimit, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));


#ifdef DEBUG
        if (result.status == EDLIB_STATUS_OK)
        {
          std::cout << "INFO, align::Aligner::doAlignment, edlib execution status = OKAY" << std::endl;
          std::cout << "INFO, align::Aligner::doAlignment, edlib execution finished, alignment length = " <<  result.alignmentLength << std::endl;
        }
        else
          std::cout << "INFO, align::Aligner::doAlignment, edlib execution status = FAILED" << std::endl;

        std::chrono::duration<double> timeAlign = skch::Time::now() - t0;
        std::cout << "INFO, align::Aligner::doAlignment, time spent= " << timeAlign.count()  << " sec" << std::endl;
#endif

        std::stringstream output;
        //Output to file
        if (result.status == EDLIB_STATUS_OK && result.alignmentLength != 0) 
        {
            uint64_t matches = 0;
            uint64_t mismatches = 0;
            uint64_t insertions = 0;
            uint64_t deletions = 0;
            uint64_t softclips = 0;
            uint64_t refAlignedLength = 0;
            uint64_t qAlignedLength = 0;

            char* cigar = alignmentToCigar(result.alignment,
                                           result.alignmentLength,
                                           refAlignedLength,
                                           qAlignedLength,
                                           matches,
                                           mismatches,
                                           insertions,
                                           deletions,
                                           softclips);

            size_t alignmentRefPos = currentRecord.rStartPos + result.startLocations[0];
            double total = refAlignedLength + (qAlignedLength - softclips);
            double identity = (double)(total - mismatches * 2 - insertions - deletions) / total;

            output << currentRecord.qId
                   << "\t" << querySize
                   << "\t" << currentRecord.qStartPos
                   << "\t" << currentRecord.qStartPos + qAlignedLength
                   << "\t" << (currentRecord.strand == skch::strnd::FWD ? "+" : "-")
                   << "\t" << refId
                   << "\t" << refSize
                   << "\t" << alignmentRefPos
                   << "\t" << alignmentRefPos + refAlignedLength
                   << "\t" << matches
                   << "\t" << std::max(refAlignedLength, qAlignedLength)
                   << "\t" << std::round(float2phred(1.0-identity))
                   << "\t" << "id:f:" << identity
                   << "\t" << "ma:i:" << matches
                   << "\t" << "mm:i:" << mismatches
                   << "\t" << "ni:i:" << insertions
                   << "\t" << "nd:i:" << deletions
                   << "\t" << "ns:i:" << softclips
                   << "\t" << "ed:i:" << result.editDistance
                   << "\t" << "al:i:" << result.alignmentLength
                   << "\t" << "se:f:" << result.editDistance * 1.0/result.alignmentLength
                   << "\t" << "cg:Z:" << cigar
                   << "\n";

            free(cigar);
        }

        edlibFreeAlignResult(result);
        delete [] queryRegionStrand;

        return output.str();
      }
  };
}


#endif
