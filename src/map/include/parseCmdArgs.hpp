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
    cmd.setIntroductoryDescription("Approximate read mapper based on Jaccard similarity");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("subject", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subject","s");

    cmd.defineOption("subjectList", "a file containing list of reference files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subjectList","sl");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("queryList", "a file containing list of query files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("queryList","ql");

    cmd.defineOption("kmer", "kmer size <= 16 [default 16 (DNA), 5 (AA)]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("pval", "p-value cutoff, used to determine window/sketch sizes [default e-03]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("pval","p");

    cmd.defineOption("window", "window size [default : computed using pvalue cutoff]\n\
P-value is not considered if a window value is provided. Lower window size implies denser sketch", 
                    ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("window","w");

    cmd.defineOption("minMatchLen", "minimum match length [default : 10000]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("minMatchLen","m");

    cmd.defineOption("perc_identity", "threshold for identity [default : 85]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("perc_identity","pi");

    cmd.defineOption("protein", "set alphabet type to proteins, default is nucleotides");
    cmd.defineOptionAlternative("protein","a");

    cmd.defineOption("filter_mode", "filter modes in mashmap: 'map', 'dot' or 'none' [default: map]\n\
'map' filters best mappings for each query sequence\n\
'dot' filters best mappings first for query and then reference sequence\n\
'none' disables filtering", 
                    ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("filter_mode", "f");

    cmd.defineOption("split", "enable split read mapping");

    cmd.defineOption("output", "output file name", ArgvParser::OptionRequired | ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");
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
   * @brief                     validate the reference and query file(s)
   * @param[in] querySequences  vector containing query file names
   * @param[in] refSequences    vector containing reference file names
   */
  template <typename VEC>
    void validateInputFiles(VEC &querySequences, VEC &refSequences)
    {
      //Open file one by one
      for(auto &e : querySequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << std::endl;
          exit(1);
        }
      }

      for(auto &e : refSequences)
      {
        std::ifstream in(e);

        if (in.fail())
        {
          std::cerr << "ERROR, skch::validateInputFiles, Could not open " << e << std::endl;
          exit(1);
        }
      }
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
    std::cout << "Match length >= " << parameters.minMatchLength << (parameters.split ? " --split":"") << std::endl;
    std::cout << "Alphabet = " << (parameters.alphabetSize == 4 ? "DNA" : "AA") << std::endl;
    std::cout << "P-value = " << parameters.p_value << std::endl;
    std::cout << "Percentage identity threshold = " << parameters.percentageIdentity << std::endl;
    std::cout << "Mapping output file = " << parameters.outFileName << std::endl;
    std::cout << "Filter mode = " << parameters.filterMode << " (1 = map, 2 = dot, 3 = none)" << std::endl;
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

    if(cmd.foundOption("protein"))
    {
      parameters.alphabetSize = 20;
    }
    else
      parameters.alphabetSize = 4;

    if(cmd.foundOption("filter_mode"))
    {
      str << cmd.optionValue("filter_mode");

      std::string filter_input;
      str >> filter_input;

      if (filter_input == "map") parameters.filterMode = filter::MAP;
      else if (filter_input == "dot") parameters.filterMode = filter::DOT;
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

    if(cmd.foundOption("split"))
    {
      parameters.split = true;
    }
    else
      parameters.split = false;


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

    if(cmd.foundOption("pval"))
    {
      str << cmd.optionValue("pval");
      str >> parameters.p_value;
      str.clear();
    }
    else
      parameters.p_value = 1e-03;

    if(cmd.foundOption("minMatchLen"))
    {
      str << cmd.optionValue("minMatchLen");
      str >> parameters.minMatchLength;
      str.clear();
    }
    else
      parameters.minMatchLength = 10000;

    if(cmd.foundOption("perc_identity"))
    {
      str << cmd.optionValue("perc_identity");
      str >> parameters.percentageIdentity;
      str.clear();
    }
    else
      parameters.percentageIdentity = 85;

    /*
     * Compute window size for sketching
     */

    if(cmd.foundOption("window"))
    {
      str << cmd.optionValue("window");
      str >> parameters.windowSize;
      str.clear();

      //Re-estimate p value
      
      int lengthQuery = parameters.split ? parameters.minMatchLength/2 : parameters.minMatchLength;

      int s = lengthQuery * 2 / parameters.windowSize; 
      parameters.p_value = skch::Stat::estimate_pvalue (s, parameters.kmerSize, parameters.alphabetSize, 
          parameters.percentageIdentity, 
          lengthQuery, parameters.referenceSize);
    }
    else
    {
      //Compute optimal window size
      parameters.windowSize = skch::Stat::recommendedWindowSize(parameters.p_value,
          parameters.kmerSize, parameters.alphabetSize,
          parameters.percentageIdentity,
          parameters.minMatchLength, parameters.referenceSize,
          parameters.split);
    }

    str << cmd.optionValue("output");
    str >> parameters.outFileName;
    str.clear();

    printCmdOptions(parameters);

    //Check if files are valid
    validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}


#endif
