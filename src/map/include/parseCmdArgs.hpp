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

    cmd.defineOption("blockLength", "keep merged mappings supported by homologies of this total length [default: segmentLength]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("blockLength", "l");

    cmd.defineOption("chainGap", "chain mappings closer than this distance in query and target, retaining mappings in best chain [default: 2*segmentLength]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("chainGap", "c");

    cmd.defineOption("numMappingsForSegment", "number of mappings to retain for each segment [default: 1]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("numMappingsForSegment", "n");

    cmd.defineOption("numMappingsForShortSeq", "number of mappings to retain for each sequence shorter than segment length [default: 1]", ArgvParser::OptionRequiresValue);
    
    cmd.defineOption("indexSave", "Filename of index file to save. Filenames w/ .tsv extension will be saved as TSV files and otherwise binary files", ArgvParser::OptionRequiresValue);
    cmd.defineOption("indexLoad", "Filename of index file to load", ArgvParser::OptionRequiresValue);


    cmd.defineOption("noSplit", "disable splitting of input sequences during mapping [enabled by default]");

    cmd.defineOption("perc_identity", "threshold for identity [default : 85]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("perc_identity","pi");

    cmd.defineOption("dropLowMapId", "drop mappings with estimated identity below --map-pct-id=\%", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("dropLowMapId", "K");

    cmd.defineOption("threads", "count of threads for parallel execution [default : 1]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("threads","t");

    cmd.defineOption("output", "output file name [default : /dev/stdout ]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("output","o");

    cmd.defineOption("kmer", "kmer size [default : 19]", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("kmer","k");

    cmd.defineOption("kmerThreshold", "ignore the top \% most-frequent kmer window [default: 0.001]", ArgvParser::OptionRequiresValue);

    cmd.defineOption("sketchSize", "Number of sketch elements", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("sketchSize","J");

    cmd.defineOption("shortenCandidateRegions", "Only compute rolling minhash score for small regions around positions where the intersection of reference and query minmers is locally maximal. Results in slighty faster runtimes at the cost of mapping placement and ANI prediction.");

    cmd.defineOption("noHgFilter", "Don't use the hypergeometric filtering and instead use the MashMap2 first pass filtering.");
    cmd.defineOption("hgFilterAniDiff", "Filter out mappings unlikely to be this ANI less than the best mapping [default: 0.0]", ArgvParser::OptionRequiresValue);
    cmd.defineOption("hgFilterConf", "Confidence value for the hypergeometric filtering [default: 0.999]", ArgvParser::OptionRequiresValue);
    cmd.defineOption("maxL1", "Maximum number of candidate regions to consider [default: unlimited]", ArgvParser::OptionRequiresValue);
    cmd.defineOption("approxMH", "");

    cmd.defineOption("skipSelf", "skip self mappings when the query and target name is the same (for all-vs-all mode)");
    cmd.defineOptionAlternative("skipSelf", "X");

    cmd.defineOption("skipPrefix", "skip mappings when the query and target have the same prefix before the last occurrence of the given character C", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("skipPrefix", "Y");

    cmd.defineOption("sparsifyMappings", "keep this fraction of mappings", ArgvParser::OptionRequiresValue);
    cmd.defineOptionAlternative("sparsifyMappings", "x");

    cmd.defineOption("noMerge", "don't merge consecutive segment-level mappings");
    cmd.defineOptionAlternative("noMerge", "M");

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
    std::cerr << "[mashmap] Reference = " << parameters.refSequences << std::endl;
    std::cerr << "[mashmap] Query = " << parameters.querySequences << std::endl;
    std::cerr << "[mashmap] Kmer size = " << parameters.kmerSize << std::endl;
    std::cerr << "[mashmap] Sketch size = " << parameters.sketchSize << std::endl;
    std::cerr << "[mashmap] Segment length = " << parameters.segLength << (parameters.split ? " (read split allowed)": " (read split disabled)") << std::endl;
    std::cerr << "[mashmap] Block length min = " << parameters.block_length << std::endl;
    std::cerr << "[mashmap] Chaining gap max = " << parameters.chain_gap << std::endl;
    std::cerr << "[mashmap] Mappings per segment = " << parameters.numMappingsForSegment << std::endl;
    std::cerr << "[mashmap] Percentage identity threshold = " << 100 * parameters.percentageIdentity << "\%" << std::endl;
    std::cerr << "[mashmap] " << (parameters.skip_self ? "Skip" : "Do not skip") << " self mappings" << std::endl;
    std::cerr << "[mashmap] " << (parameters.stage2_full_scan ? "Full scan" : "Short scan") << std::endl;
    if (parameters.stage1_topANI_filter) 
      std::cerr << "[mashmap] " << "Hypergeometric filter w/ delta = " << parameters.ANIDiff << " and confidence " << parameters.ANIDiffConf << std::endl;
    else
      std::cerr << "[mashmap] " <<  "No hypergeometric filter" << std::endl;
    if (parameters.maxL1) 
      std::cerr << "[mashmap] " <<  "At most " << parameters.maxL1 << " L1 candidate regions per segment" << std::endl;
    else
      std::cerr << "[mashmap] " <<  "Unlimited L1 candidate regions per segment" << std::endl;

    std::cerr << "[mashmap] Mapping output file = " << parameters.outFileName << std::endl;
    std::cerr << "[mashmap] Filter mode = " << parameters.filterMode << " (1 = map, 2 = one-to-one, 3 = none)" << std::endl;
    std::cerr << "[mashmap] Execution threads  = " << parameters.threads << std::endl;
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
    else if (cmd.foundOption("queryList"))
    {
      std::string listFile;

      str << cmd.optionValue("queryList");
      str >> listFile;

      parseFileList(listFile, parameters.querySequences);
    } 
    else
    {
      parameters.skip_self = true;
      parameters.querySequences = parameters.refSequences;
    }
    str.clear();

    if (cmd.foundOption("skipSelf"))
    {
        parameters.skip_self = true;
    } else {
        parameters.skip_self = false;
    }

    parameters.approximateMinHash = cmd.foundOption("approxMH");

    if (cmd.foundOption("skipPrefix"))
    {
      str << cmd.optionValue("skipPrefix");
      str >> parameters.prefix_delim;
      parameters.skip_prefix = true;
    } else {
      parameters.skip_prefix = false;
      parameters.prefix_delim = '\0';
    }

    str.clear();

    if (cmd.foundOption("indexSave")) {
        str << cmd.optionValue("indexSave");
        str >> parameters.saveIndexFilename;
    } else {
        parameters.saveIndexFilename = "";
    }
    if (cmd.foundOption("indexLoad")) {
        str << cmd.optionValue("indexLoad");
        str >> parameters.loadIndexFilename;
    } else {
        parameters.loadIndexFilename = "";
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

    if(cmd.foundOption("noMerge"))
    {
      parameters.mergeMappings = false;
    }
    else
      parameters.mergeMappings = true;


    //Parse algorithm parameters
    if(cmd.foundOption("kmer"))
    {
      str << cmd.optionValue("kmer");
      str >> parameters.kmerSize;
      str.clear();
    }
    else
    {
      parameters.kmerSize = 19;
    }
    str.clear();

    if(cmd.foundOption("segLength"))
    {
      str << cmd.optionValue("segLength");
      str >> parameters.segLength;
      str.clear();

      if(parameters.segLength < 100)
      {
        std::cerr << "ERROR, skch::parseandSave, minimum segment length is required to be >= 100 bp.\n\
          This is because Mashmap is not designed for computing short local alignments.\n" << std::endl;
        exit(1);
      }
    }
    else
      parameters.segLength = 5000;

    if(cmd.foundOption("blockLength"))
    {
      str << cmd.optionValue("blockLength");
      str >> parameters.block_length;
      if (parameters.block_length < 0) {
          std::cerr << "[mashmap] ERROR, skch::parseandSave, min block length has to be a float value greater than or equal to 0." << std::endl;
          exit(1);
      }
    } else {
        // n.b. we map-merge across gaps up to 3x segment length
        // and then filter for things that are at least block_length long
        parameters.block_length = parameters.segLength;
    }
    str.clear();

    if (cmd.foundOption("chainGap")) {
      int64_t l;
      str << cmd.optionValue("chainGap");
      str >> l;
      if (l < 0) {
          std::cerr << "[mashmap] ERROR, skch::parseandSave, chain gap has to be a float value greater than or equal to 0." << std::endl;
          exit(1);
      }
      parameters.chain_gap = l;
    } else {
      parameters.chain_gap = 2*parameters.segLength;
    }
    str.clear();

    if (cmd.foundOption("dropLowMapId")) {
        parameters.keep_low_pct_id = false;
    } else {
        parameters.keep_low_pct_id = true;
    }

    if (cmd.foundOption("kmerThreshold")) {
        str << cmd.foundOption("kmerThreshold");
        str >> parameters.kmer_pct_threshold;
    } else {
        parameters.kmer_pct_threshold = 0.001; // in percent! so we keep 99.999% of kmers
    }
    str.clear();

    if (cmd.foundOption("numMappingsForSegment")) {
      uint32_t n;
      str << cmd.optionValue("numMappingsForSegment");
      str >> n;
      if (n > 0) {
          parameters.numMappingsForSegment = n;
      } else {
          std::cerr << "[mashmap] ERROR, skch::parseandSave, the number of mappings to retain for each segment has to be greater than 0." << std::endl;
          exit(1);
      }
    } else {
      parameters.numMappingsForSegment = 1;
    }
    str.clear();

    if (cmd.foundOption("numMappingsForShortSeq")) {
      uint32_t n;
      str << cmd.optionValue("numMappingsForShortSeq");
      str >> n;
      if (n > 0) {
          parameters.numMappingsForShortSequence = n;
      } else {
          std::cerr << "[mashmap] ERROR, skch::parseandSave, the number of mappings to retain for each sequence shorter than segment length has to be grater than 0." << std::endl;
          exit(1);
      }
    } else {
        parameters.numMappingsForShortSequence = 1;
    }
    str.clear();

    if(cmd.foundOption("perc_identity"))
    {
      str << cmd.optionValue("perc_identity");
      str >> parameters.percentageIdentity;
      str.clear();

      if(parameters.percentageIdentity < 50)
      {
        std::cerr << "ERROR, skch::parseandSave, minimum nucleotide identity requirement should be >= 50\%\n" << std::endl;
        exit(1);
      }
      parameters.percentageIdentity /= 100.0;
    }
    else
      parameters.percentageIdentity = 0.85;
    str.clear();

    parameters.stage1_topANI_filter = !cmd.foundOption("noHgFilter"); 

    if (cmd.foundOption("hgFilterAniDiff")) {
      str << cmd.optionValue("hgFilterAniDiff");
      str >> parameters.ANIDiff;
      if (parameters.ANIDiff < 0 ||  parameters.ANIDiff > 100) {
        std::cerr << "ERROR, skch::parseandSave, ANI difference must be between 0 and 100" << std::endl;
        exit(1);
      }
    } else {
      parameters.ANIDiff = skch::fixed::ANIDiff;
    }
    str.clear();

    if (cmd.foundOption("hgFilterConf")) {
      str << cmd.optionValue("hgFilterConf");
      str >> parameters.ANIDiffConf;
      if (parameters.ANIDiffConf < 0 ||  parameters.ANIDiffConf > 1) {
        std::cerr << "ERROR, skch::parseandSave, hypergeometric confidence must be between 0 and 1" << std::endl;
        exit(1);
      }
    } else {
      parameters.ANIDiffConf = skch::fixed::ANIDiffConf;
    }
    str.clear();

    if (cmd.foundOption("maxL1")) {
      str << cmd.optionValue("maxL1");
      str >> parameters.maxL1;
    }
    else 
      parameters.maxL1 = 0;
    str.clear();

    parameters.stage2_full_scan = !cmd.foundOption("shortenCandidateRegions");

    if (cmd.foundOption("sparsifyMappings")) {
      double frac;
      str << cmd.optionValue("sparsifyMappings");
      str >> frac;
      if (frac == 1) {
          // overflows
          parameters.sparsity_hash_threshold = std::numeric_limits<uint64_t>::max();
      } else {
          parameters.sparsity_hash_threshold
              = frac * std::numeric_limits<uint64_t>::max();
      }
    } else {
        parameters.sparsity_hash_threshold
            = std::numeric_limits<uint64_t>::max();
    }
    str.clear();

    if(cmd.foundOption("threads"))
    {
      str << cmd.optionValue("threads");
      str >> parameters.threads;
      str.clear();
    }
    else
      parameters.threads = 1;


    //Compute the best sketch size
    if(cmd.foundOption("sketchSize")) 
    {
      str << cmd.optionValue("sketchSize");
      str >> parameters.sketchSize;
      str.clear();
    } else {
      if (parameters.percentageIdentity >= 0.95) {
        parameters.sketchSize = parameters.segLength / 25;
      }
      else if (parameters.percentageIdentity >= 0.85) {
        parameters.sketchSize = parameters.segLength / 16;
      }
      else if (parameters.percentageIdentity >= 0.75) {
        parameters.sketchSize = parameters.segLength / 13;
      } 
      else {
        parameters.sketchSize = parameters.segLength / 10;
      }
    }

    if(cmd.foundOption("output"))
    {
      str << cmd.optionValue("output");
      str >> parameters.outFileName;
      str.clear();
    }
    else
      parameters.outFileName = "/dev/stdout";


    printCmdOptions(parameters);

    //Check if files are valid
    validateInputFiles(parameters.querySequences, parameters.refSequences);
  }
  

}


#endif
