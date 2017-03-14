/**
 * @file    winmin_parseCmdArgs.hpp
 * @brief   Functionality related to command line parsing for indexing and mapping
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef WINMIN_PARSE_CMD_HPP 
#define WINMIN_PARSE_CMD_HPP

#include <iostream>
#include <string>
#include <fstream>

//Own includes
#include "winnowed_minhash/include/parameters.hpp"

//External includes
#include "argvparser.hpp"

namespace winmin
{
  /**
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd  command line parser object
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd)
  {
    cmd.setIntroductoryDescription("Approximate read mapper based on Jaccard similarity");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("subject", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOptionAlternative("subject","s");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("kmer", "kmer size <= 16 [default 16]", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("window", "window size", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);
    cmd.defineOptionAlternative("window","w");
  }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(Parameters &parameters)
  {
    std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "Reference = " << parameters.refSequence << std::endl;
    std::cout << "Query = " << parameters.querySequence << std::endl;
    std::cout << "Kmer size = " << parameters.kmerSize << std::endl;
    std::cout << "Window size = " << parameters.windowSize << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  sketch parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      CommandLineProcessing::ArgvParser &cmd, 
      Parameters &parameters)
  {
    int result = cmd.parse(argc, argv);

    //Make sure we get the right command line args
    if (result != ArgvParser::NoParserError)
    {
      std::cout << cmd.parseErrorDescription(result) << "\n";
      exit(1);
    }

    std::stringstream str;

    //Parse reference files
    if(cmd.foundOption("subject"))
    {
      std::string ref;

      str << cmd.optionValue("subject");
      str >> ref;

      parameters.refSequence = ref;
    }

    str.clear();

    //Parse query files
    if(cmd.foundOption("query"))
    {
      std::string query;

      str << cmd.optionValue("query");
      str >> query;

      parameters.querySequence = query;
    }

    str.clear();

    //Parse algorithm parameters
    if(cmd.foundOption("kmer"))
    {
      str << cmd.optionValue("kmer");
      str >> parameters.kmerSize;
      str.clear();
    }
    else
    {
      parameters.kmerSize = 16;
    }

    /*
     * Compute window size for sketching
     */

    if(cmd.foundOption("window"))
    {
      str << cmd.optionValue("window");
      str >> parameters.windowSize;
      str.clear();
    }

    printCmdOptions(parameters);
  }
}


#endif
