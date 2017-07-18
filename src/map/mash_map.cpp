/**
 * @file    mash_map.cpp
 * @ingroup src
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>
#include <omp.h>

//Own includes
#include "map/include/map_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/parseCmdArgs.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"

//External includes
#include "common/argvparser.hpp"

int main(int argc, char** argv)
{
  CommandLineProcessing::ArgvParser cmd;

  //Setup command line options
  skch::initCmdParser(cmd);

  //Parse command line arguements   
  skch::Parameters parameters;        //sketching and mapping parameters

  skch::parseandSave(argc, argv, cmd, parameters);   

  omp_set_num_threads(parameters.threads);

  auto t0 = skch::Time::now();

  //Build the sketch for reference
  skch::Sketch referSketch(parameters);

  std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
  std::cout << "INFO, skch::main, Time spent sketching the reference : " << timeRefSketch.count() << " sec" << std::endl;

  //Map the sequences in query file
  t0 = skch::Time::now();

  skch::Map mapper = skch::Map(parameters, referSketch);

  std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
  std::cout << "INFO, skch::main, Time spent mapping the query : " << timeMapQuery.count() << " sec" << std::endl;

  std::cout << "INFO, skch::main, mapping results saved in : " << parameters.outFileName << std::endl;

}
