/**
 * @file    align.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <sstream>
#include <fstream>
#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>
#include <cstdio>

//Own includes
#include "map/include/base_types.hpp"
#include "align/include/align_parameters.hpp"
#include "align/include/parseCmdArgs.hpp"
#include "align/include/computeAlignments.hpp"

//External includes
#include "common/argvparser.hpp"

int main(int argc, char** argv)
{
  /*
   * Make sure env variable MALLOC_ARENA_MAX is unset 
   * for efficient multi-thread execution
   */
  unsetenv((char *)"MALLOC_ARENA_MAX");

  CommandLineProcessing::ArgvParser cmd;

  //Setup command line options
  align::initCmdParser(cmd);

  //Parse command line arguements   
  align::Parameters parameters;        //sketching and mapping parameters

  align::parseandSave(argc, argv, cmd, parameters);   

  auto t0 = skch::Time::now();

  align::Aligner alignObj(parameters);

  std::chrono::duration<double> timeRefRead = skch::Time::now() - t0;
  std::cout << "INFO, align::main, Time spent read the reference sequences: " << timeRefRead.count() << " sec" << std::endl;

  //Compute the alignments
  alignObj.compute();

  std::chrono::duration<double> timeAlign = skch::Time::now() - t0;
  std::cout << "INFO, align::main, Time spent computing the aligment: " << timeAlign.count() << " sec" << std::endl;

  std::cout << "INFO, align::main, alignment results saved in: " << parameters.samOutputFile << std::endl;
}
