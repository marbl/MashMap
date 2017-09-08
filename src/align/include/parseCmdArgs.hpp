/**
 * @file    parseCmdArgs.hpp
 * @brief   Functionality related to aligning mashmap mappings
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PARSE_CMD_ALIGN_HPP 
#define PARSE_CMD_ALIGN_HPP

#include <iostream>
#include <string>
#include <fstream>

//Own includes
#include "align/include/align_parameters.hpp"
#include "map/include/parseCmdArgs.hpp"

//External includes
#include "common/argvparser.hpp"

namespace align
{
  /**
   * @brief           Initialize the command line argument parser 
   * @param[out] cmd  command line parser object
   */
  void initCmdParser(CommandLineProcessing::ArgvParser &cmd)
  {
    cmd.setIntroductoryDescription("-----------------\n\
Post process mashmap output to compute alignments for obtaining SAM output.\n\
Provide same reference and query files that were used for obtaining mashmap mapping boundaries.\n\
-----------------\n\
Example usage: \n\
$ mashmap-align -s ref.fa -q seq.fq --mappingFile mashmap.out [OPTIONS]");

    cmd.setHelpOption("h", "help", "Print this help page");

    cmd.defineOption("subject", "an input reference file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subject","s");

    cmd.defineOption("subjectList", "a file containing list of reference files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("subjectList","sl");

    cmd.defineOption("query", "an input query file (fasta/fastq)[.gz]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("query","q");

    cmd.defineOption("queryList", "a file containing list of query files, one per line", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("queryList","ql");

    cmd.defineOption("mappingFile", "mashmap file containing mapping information", ArgvParser::OptionRequiresValue | ArgvParser::OptionRequired);

    cmd.defineOption("threads", "count of threads for parallel execution [default : 1]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("threads","t");

    cmd.defineOption("output", "output file name [default : mashmap.out.sam]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");
  }

  /**
   * @brief                   Print the parsed cmd line options
   * @param[in]  parameters   parameters parsed from command line
   */
  void printCmdOptions(align::Parameters &parameters)
  {
    std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
    std::cout << "Reference = " << parameters.refSequences << std::endl;
    std::cout << "Query = " << parameters.querySequences << std::endl;
    std::cout << "Mapping file = " << parameters.mashmapPafFile << std::endl;
    std::cout << "Alignment output file = " << parameters.samOutputFile << std::endl;
    std::cout << ">>>>>>>>>>>>>>>>>>" << std::endl;
  }

  /**
   * @brief                   Parse the cmd line options
   * @param[in]   cmd
   * @param[out]  parameters  alignment parameters are saved here
   */
  void parseandSave(int argc, char** argv, 
      CommandLineProcessing::ArgvParser &cmd, 
      align::Parameters &parameters)
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
      std::cerr << "ERROR, align::parseandSave, Provide reference file(s)\n\
        This input should be same as used for generating mashmap mapping output" << std::endl;
      exit(1);
    }
    else if (!cmd.foundOption("query") && !cmd.foundOption("queryList"))
    {
      std::cerr << "ERROR, align::parseandSave, Provide query file(s)\n\
        This input should be same as used for generating mashmap mapping output" << std::endl;
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

      skch::parseFileList(listFile, parameters.refSequences);
    }

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

      skch::parseFileList(listFile, parameters.querySequences);
    }

    str.clear();

    str << cmd.optionValue("mappingFile");
    str >> parameters.mashmapPafFile;

    str.clear();

    if(cmd.foundOption("threads"))
    {
      str << cmd.optionValue("threads");
      str >> parameters.threads;
    }
    else
      parameters.threads = 1;

    str.clear();

    if(cmd.foundOption("output"))
    {
      str << cmd.optionValue("output");
      str >> parameters.samOutputFile;
      str.clear();
    }
    else
      parameters.samOutputFile = "mashmap.out.sam";

    str.clear();

    printCmdOptions(parameters);

    //Check if files are valid
    skch::validateInputFile(parameters.mashmapPafFile);
    skch::validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
}

#endif
