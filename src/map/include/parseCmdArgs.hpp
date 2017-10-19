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
$ mashmap -s ref.fa -q seq.fq [OPTIONS]\n\
$ mashmap --sl reference_files_list.txt -q seq.fq [OPTIONS]");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("subject", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subject","s");

    cmd.defineOption("subjectList", "a file containing list of reference files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subjectList","sl");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("queryList", "a file containing list of query files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("queryList","ql");

    cmd.defineOption("segLength", "mapping segment length [default : 5,000]\n\
sequences shorter than segment length will be ignored", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("segLength","m");

    cmd.defineOption("noSplit", "disable splitting of input sequences during mapping [enabled by default]");

    cmd.defineOption("perc_identity", "threshold for identity [default : 85]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("perc_identity","pi");

    cmd.defineOption("threads", "count of threads for parallel execution [default : 1]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("threads","t");

    cmd.defineOption("output", "output file name [default : mashmap.out]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");

    cmd.defineOption("kmer", "kmer size <= 16 [default : 16]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("filter_mode", "filter modes in mashmap: 'map', 'one-to-one' or 'none' [default: map]\n\
'map' computes best mappings for each query sequence\n\
'one-to-one' computes best mappings for query as well as reference sequence\n\
'none' disables filtering", 
                    ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("filter_mode", "f");
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

      for(auto &e : refSequences)
        validateInputFile(e);
    }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(skch::Parameters &parameters)
  {
    std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "Reference = " << parameters.refSequences << std::endl;
    std::cout << "Query = " << parameters.querySequences << std::endl;
    std::cout << "Kmer size = " << parameters.kmerSize << std::endl;
    std::cout << "Window size = " << parameters.windowSize << std::endl;
    std::cout << "Segment length = " << parameters.segLength << (parameters.split ? " (read split allowed)": " (read split disabled)") << std::endl;
    std::cout << "Alphabet = " << (parameters.alphabetSize == 4 ? "DNA" : "AA") << std::endl;
    std::cout << "Percentage identity threshold = " << parameters.percentageIdentity << "\%" << std::endl;
    std::cout << "Mapping output file = " << parameters.outFileName << std::endl;
    std::cout << "Filter mode = " << parameters.filterMode << " (1 = map, 2 = one-to-one, 3 = none)" << std::endl;
    std::cout << "Execution threads  = " << parameters.threads << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
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
    else if (!cmd.foundOption("subject") && !cmd.foundOption("subjectList"))
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
    if(cmd.foundOption("subject"))
    {
      std::string ref;

      str << cmd.optionValue("subject");
      str >> ref;

      parameters.refSequences.push_back(ref);
    }
    else //list of files
    {
      std::string listFile;

      str << cmd.optionValue("subjectList");
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

      if(parameters.segLength < 1000)
      {
        std::cerr << "ERROR, skch::parseandSave, minimum local alignment length requirement should be >= 1000 bp.\n\
          This is because Mashmap is not efficient for computing short alignments.\n" << std::endl;
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

    /*
     * Compute window size for sketching
     */

    //Compute optimal window size
    parameters.windowSize = skch::Stat::recommendedWindowSize(skch::fixed::pval_cutoff,
        parameters.kmerSize, parameters.alphabetSize,
        parameters.percentageIdentity,
        parameters.segLength, parameters.referenceSize);

    if(cmd.foundOption("output"))
    {
      str << cmd.optionValue("output");
      str >> parameters.outFileName;
      str.clear();
    }
    else
      parameters.outFileName = "mashmap.out";


    printCmdOptions(parameters);

    //Check if files are valid
    validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}


#endif
