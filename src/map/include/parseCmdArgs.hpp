/**
 * @file    parseCmdArgs.hpp
 * @brief   Functionality related to command line parsing for indexing and mapping
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PARSE_CMD_HPP 
#define PARSE_CMD_HPP

#include <iostream>
#include <string>
#include <fstream>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/map_stats.hpp"
#include "map/include/commonFunc.hpp"

//External includes
#include "common/argvparser.hpp"

namespace skch
{

  /**
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd  command line parser object
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd)
  {
    cmd.setIntroductoryDescription("-----------------\n\
Mashmap is an approximate long read or contig mapper based on Jaccard similarity\n\
-----------------\n\
Example usage: \n\
$ mashmap -r ref.fa -q seq.fq [OPTIONS]\n\
$ mashmap --rl reference_files_list.txt -q seq.fq [OPTIONS]");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("ref", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("ref","r");

    cmd.defineOption("refList", "a file containing list of reference files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("refList","rl");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("queryList", "a file containing list of query files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("queryList","ql");

    cmd.defineOption("segLength", "mapping segment length [default : 5,000]\n\
sequences shorter than segment length will be ignored", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("segLength","s");

    cmd.defineOption("noSplit", "disable splitting of input sequences during mapping [enabled by default]");

    cmd.defineOption("perc_identity", "threshold for identity [default : 85]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("perc_identity","pi");

    cmd.defineOption("threads", "count of threads for parallel execution [default : 1]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("threads","t");

    cmd.defineOption("output", "output file name [default : mashmap.out]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");

    cmd.defineOption("kmer", "kmer size <= 16 [default : 16]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");
    
    cmd.defineOption("sketchSize", "size of sketch [default : 25]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("sketchSize","S");

    cmd.defineOption("filter_mode", "filter modes in mashmap: 'map', 'one-to-one' or 'none' [default: map]\n\
'map' computes best mappings for each query sequence\n\
'one-to-one' computes best mappings for query as well as reference sequence\n\
'none' disables filtering", 
                    ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("filter_mode", "f");

    cmd.defineOption("secondaries", "number of secondary mappings in 'map' filter_mode [default : 0]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("secondaries", "n");
  }

  /**
   * @brief                   Parse the file which has list of reference or query files
   * @param[in]   fileToRead  File containing list of ref/query files 
   * @param[out]  fileList    List of files will be saved in this vector   
   */
  template <typename VEC>
    void parseFileList(std::string &fileToRead, VEC &fileList)
    {
      std::string line;

      std::ifstream in(fileToRead);

      if (in.fail())
      {
        std::cerr << "ERROR, skch::parseFileList, Could not open " << fileToRead << std::endl;
        exit(1);
      }

      while (std::getline(in, line))
      {
        fileList.push_back(line);
      }
    }

  /**
   * @brief                     validate single input file
   * @param[in] fileName        file name
   */
  void validateInputFile(std::string &fileName)
  {
    //Open file one by one
    std::ifstream in(fileName);

    if (in.fail())
    {
      std::cerr << "ERROR, skch::validateInputFile, Could not open " << fileName << std::endl;
      exit(1);
    }
  }

  bool ends_with_string(std::string const& str, std::string const& what) {
      return what.size() <= str.size()
      && str.find(what, str.size() - what.size()) != str.npos;
  }

  bool checkIndexFileExists(std::string fastaName, const char* suffix){
      std::string indexFileName(fastaName);
      indexFileName = indexFileName + suffix;
      std::ifstream in(indexFileName);
      return !in.fail();
  }

  /**
 * @brief                     check if the indexes for the input file exist (.fai for fasta, .fai + .gzi for bgzipped fasta)
 * @param[in] fileName        file name
 */
  void validateInputFileIndexes(std::string &fileName)
  {
      bool indexesExist = true;

      if (ends_with_string(fileName, ".fa") || ends_with_string(fileName, ".fasta") || ends_with_string(fileName, ".fna")) {
          indexesExist &= checkIndexFileExists(fileName, ".fai");
      }else if (ends_with_string(fileName, ".fa.gz") || ends_with_string(fileName, ".fasta.gz") || ends_with_string(fileName, ".fna.gz")) {
          indexesExist &= checkIndexFileExists(fileName, ".fai");
          indexesExist &= checkIndexFileExists(fileName, ".gzi");
      } else {
          std::cerr << "\tUnknown type for file: "<< fileName << std::endl;
          exit(1);
      }

      if (!indexesExist)
      {
          std::cerr << "ERROR, skch::validateInputFileIndexes, missing index(es) for the file "<< fileName << std::endl;
          if (ends_with_string(fileName, ".fa") || ends_with_string(fileName, ".fasta") || ends_with_string(fileName, ".fna")) {
              std::cerr << "\tWe recommend to build the indexes on BGZIP files by executing:" << std::endl;
              std::cerr << "\tbgzip -@ number_of_threads "<< fileName << std::endl;
              std::cerr << "\tsamtools faidx "<< fileName << ".gz" << std::endl;
          } else {
              std::cerr << "\tIf the compressed file is in BGZIP, execute:" << std::endl;
              std::cerr << "\tsamtools faidx "<< fileName << std::endl;
          }

          exit(1);
      }
  }

  /**
   * @brief                     validate the reference and query file(s)
   * @param[in] querySequences  vector containing query file names
   * @param[in] refSequences    vector containing reference file names
   */
  template <typename VEC>
    void validateInputFiles(VEC &querySequences, VEC &refSequences)
    {
      //validate files one by one
      for(auto &e : querySequences)
        validateInputFile(e);

      for(auto &e : refSequences) {
          validateInputFile(e);
          validateInputFileIndexes(e);
      }
    }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(skch::Parameters &parameters)
  {
    std::cerr << "[mashmap::map] Reference = " << parameters.refSequences << std::endl;
    std::cerr << "[mashmap::map] Query = " << parameters.querySequences << std::endl;
    std::cerr << "[mashmap::map] Kmer size = " << parameters.kmerSize << std::endl;
    std::cerr << "[mashmap::map] Sketch size = " << parameters.sketchSize << std::endl;
    std::cerr << "[mashmap::map] Segment length = " << parameters.segLength << (parameters.split ? " (read split allowed)": " (read split disabled)") << std::endl;
    std::cerr << "[mashmap::map] Block length min = " << parameters.block_length << std::endl;
    std::cerr << "[mashmap::map] Chaining gap max = " << parameters.chain_gap << std::endl;
    //std::cerr << "[mashmap::map] Alphabet = " << (parameters.alphabetSize == 4 ? "DNA" : "AA") << std::endl;
    std::cerr << "[mashmap::map] Percentage identity threshold = " << 100 * parameters.percentageIdentity << "\%" << std::endl;
    std::cerr << "[mashmap::map] " << (parameters.skip_self ? "Skip" : "Do not skip") << " self mappings" << std::endl;
    std::cerr << "[mashmap::map] Mapping output file = " << parameters.outFileName << std::endl;
    std::cerr << "[mashmap::map] Filter mode = " << parameters.filterMode << " (1 = map, 2 = one-to-one, 3 = none)" << std::endl;
    std::cerr << "[mashmap::map] Execution threads  = " << parameters.threads << std::endl;
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      CommandLineProcessing::ArgvParser &cmd, 
      skch::Parameters &parameters)
  {
    int result = cmd.parse(argc, argv);

    //Make sure we get the right command line args
    if (result != ArgvParser::NoParserError)
    {
      std::cerr << cmd.parseErrorDescription(result) << std::endl;
      exit(1);
    }
    else if (!cmd.foundOption("ref") && !cmd.foundOption("refList"))
    {
      std::cerr << "ERROR, skch::parseandSave, Provide reference file(s)" << std::endl;
      exit(1);
    }
    else if (!cmd.foundOption("query") && !cmd.foundOption("queryList"))
    {
      std::cerr << "ERROR, skch::parseandSave, Provide query file(s)" << std::endl;
      exit(1);
    }

    std::stringstream str;

    //Parse reference files
    if(cmd.foundOption("ref"))
    {
      std::string ref;

      str << cmd.optionValue("ref");
      str >> ref;

      parameters.refSequences.push_back(ref);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("refList");
      str >> listFile;

      parseFileList(listFile, parameters.refSequences);
    }

    //Size of reference
    parameters.referenceSize = skch::CommonFunc::getReferenceSize(parameters.refSequences); 
    str.clear();

    //Parse query files
    if(cmd.foundOption("query"))
    {
      std::string query;

      str << cmd.optionValue("query");
      str >> query;

      parameters.querySequences.push_back(query);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("queryList");
      str >> listFile;

      parseFileList(listFile, parameters.querySequences);
    }

    str.clear();

    parameters.alphabetSize = 4;
    //Do not expose the option to set protein alphabet in mashmap
    //parameters.alphabetSize = 20;

    if(cmd.foundOption("filter_mode"))
    {
      str << cmd.optionValue("filter_mode");

      std::string filter_input;
      str >> filter_input;

      if (filter_input == "map") parameters.filterMode = filter::MAP;
      else if (filter_input == "one-to-one") parameters.filterMode = filter::ONETOONE;
      else if (filter_input == "none") parameters.filterMode = filter::NONE;
      else 
      {
        std::cerr << "ERROR, skch::parseandSave, Invalid option given for filter_mode" << std::endl;
        exit(1);
      };

      str.clear();
    }
    else
      parameters.filterMode = filter::MAP;

    if(cmd.foundOption("noSplit"))
    {
      parameters.split = false;
    }
    else
      parameters.split = true;

    if(cmd.foundOption("sketchSize"))
    {
      str << cmd.optionValue("sketchSize");
      str >> parameters.sketchSize;
      str.clear();
    }

    //Parse algorithm parameters
    if(cmd.foundOption("kmer"))
    {
      str << cmd.optionValue("kmer");
      str >> parameters.kmerSize;
      str.clear();
    }
    else
    {
      if(parameters.alphabetSize == 4)
        parameters.kmerSize = 16;
      else
        parameters.kmerSize = 5;
    }

    if(cmd.foundOption("segLength"))
    {
      str << cmd.optionValue("segLength");
      str >> parameters.segLength;
      str.clear();

      if(parameters.segLength < 500)
      {
        std::cerr << "ERROR, skch::parseandSave, minimum segment length is required to be >= 500 bp.\n\
          This is because Mashmap is not designed for computing short local alignments.\n" << std::endl;
        exit(1);
      }
    }
    else
      parameters.segLength = 5000;

    if(cmd.foundOption("perc_identity"))
    {
      str << cmd.optionValue("perc_identity");
      str >> parameters.percentageIdentity;
      str.clear();

      if(parameters.percentageIdentity < 70)
      {
        std::cerr << "ERROR, skch::parseandSave, minimum nucleotide identity requirement should be >= 70\%\n" << std::endl;
        exit(1);
      }
    }
    else
      parameters.percentageIdentity = 85;

    if(cmd.foundOption("threads"))
    {
      str << cmd.optionValue("threads");
      str >> parameters.threads;
      str.clear();
    }
    else
      parameters.threads = 1;


    if(cmd.foundOption("output"))
    {
      str << cmd.optionValue("output");
      str >> parameters.outFileName;
      str.clear();
    }
    else
      parameters.outFileName = "mashmap.out";

    if(cmd.foundOption("secondaries"))
    {
      str << cmd.optionValue("secondaries");
      str >> parameters.numMappingsForSegment;
      str.clear();
    }
    else
    {
      parameters.numMappingsForSegment = 1;
    }

    printCmdOptions(parameters);

    //Check if files are valid
    validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}


#endif
