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

#include "map/include/map_parameters.hpp"
#include "map/include/base_types.hpp"
#include "map/include/winSketch.hpp"
#include "map/include/computeMap.hpp"
#include "map/include/parseCmdArgs.hpp"

#include "interface/parse_args.hpp"


//External includes
#include "common/args.hxx"
#include "common/ALeS.hpp"

int main(int argc, char** argv) {
    /*
     * Make sure env variable MALLOC_ARENA_MAX is unset 
     * for efficient multi-thread execution
     */
    unsetenv((char *)"MALLOC_ARENA_MAX");

    // get our parameters from the command line
    skch::Parameters map_parameters;
    yeet::Parameters yeet_parameters;
    yeet::parse_args(argc, argv, map_parameters, yeet_parameters);

    //parameters.refSequences.push_back(ref);

    //skch::parseandSave(argc, argv, cmd, parameters);
    skch::printCmdOptions(map_parameters);

    auto t0 = skch::Time::now();

    //Build the sketch for reference
    skch::Sketch referSketch(map_parameters);

    std::chrono::duration<double> timeRefSketch = skch::Time::now() - t0;
    std::cerr << "[mashmap::map] time spent computing the reference index: " << timeRefSketch.count() << " sec" << std::endl;

    //Map the sequences in query file
    t0 = skch::Time::now();

    skch::Map mapper = skch::Map(map_parameters, referSketch);

    std::chrono::duration<double> timeMapQuery = skch::Time::now() - t0;
    std::cerr << "[mashmap::map] time spent mapping the query: " << timeMapQuery.count() << " sec" << std::endl;
    std::cerr << "[mashmap::map] mapping results saved in: " << map_parameters.outFileName << std::endl;

    return 0;
}
