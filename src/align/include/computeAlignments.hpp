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
#include <zlib.h>  
#include <cassert>

//Own includes
#include "align/include/align_types.hpp"
#include "align/include/align_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/kseq.h"
#include "common/edlib.h"

KSEQ_INIT(gzFile, gzread)

namespace align
{
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

      /**
       * @brief                 parse query sequences and mashmap mappings
       *                        to compute sequence alignments
       */
      void computeAlignments()
      {
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

          std::ofstream outstrm(param.samOutputFile);

          while ((len = kseq_read(seq)) >= 0) 
          {
             std::string qSeqId = seq->name.s;

             skch::CommonFunc::makeUpperCase(seq->seq.s, len);
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
               this->doAlignment(currentRecord, mappingRecordLine, qSequence, outstrm);
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
                 this->doAlignment(currentRecord, mappingRecordLine, qSequence, outstrm);
               }
             }
          }

          mappingListStream.close();

          kseq_destroy(seq);  
          gzclose(fp); //close the file handler 
          fclose(file);
        }
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
          currentRecord.strand =  tokens[4] == "+" ? skch::strnd::FWD : skch::strnd::REV;;
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
      void doAlignment(MappingBoundaryRow &currentRecord, const std::string &mappingRecordLine, const std::string &qSequence, std::ofstream &outstrm)
      {

#ifdef DEBUG
        std::cout << "INFO, align::Aligner::doAlignment, aligning mashmap record: " << mappingRecordLine << std::endl;
#endif

        //Define reference substring for this mapping
        const std::string &refId = currentRecord.refId;
        const char* refRegion = this->refSequences[refId].c_str();
        refRegion += currentRecord.rStartPos;
        skch::offset_t refLen = currentRecord.rEndPos - currentRecord.rStartPos + 1;

        assert(refLen <= this->refSequences[refId].size());

        //Define query substring for this mapping
        const char* queryRegion = qSequence.c_str();  //initially point to beginning
        skch::offset_t queryLen = currentRecord.qEndPos - currentRecord.qStartPos + 1;
        queryRegion += currentRecord.qStartPos;

        char* queryRegionStrand = new char[queryLen];

        if(currentRecord.strand == skch::strnd::FWD) 
          strncpy(queryRegionStrand, queryRegion, queryLen);    //Copy the same string
        else
          skch::CommonFunc::reverseComplement(queryRegion, queryRegionStrand, queryLen); //Reverse complement

        assert(queryLen <= qSequence.size());

        //Compute alignment

        EdlibAlignResult result = edlibAlign(queryRegionStrand, queryLen, refRegion, refLen, 
            edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));

#ifdef DEBUG
        if (result.status == EDLIB_STATUS_OK)
          std::cout << "INFO, align::Aligner::doAlignment, edlib execution status = OKAY" << std::endl;
        else
          std::cout << "INFO, align::Aligner::doAlignment, edlib execution status = FAILED" << std::endl;
#endif

        //Output to file
        if (result.status == EDLIB_STATUS_OK) 
        {
          char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);

          outstrm << mappingRecordLine 
            << " " << result.editDistance * 1.0/result.alignmentLength
            << " " << cigar 
            << "\n";

          free(cigar);
        }

        edlibFreeAlignResult(result);
        delete [] queryRegionStrand;
      }
  };
}


#endif
