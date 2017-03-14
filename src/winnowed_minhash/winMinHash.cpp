/**
 * @file    winMinHash.cpp
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#include <iostream>
#include <ctime>
#include <chrono>
#include <functional>

//Own includes
#include "winnowed_minhash/include/parameters.hpp"
#include "winnowed_minhash/include/parseCmdArgs.hpp"
#include "winnowed_minhash/include/computeJaccard.hpp"

//External includes
#include "argvparser.hpp"

int main(int argc, char** argv)
{
  CommandLineProcessing::ArgvParser cmd;

  //Setup command line options
  winmin::initCmdParser(cmd);

  //Parse command line arguements   
  winmin::Parameters parameters;        //sketching and mapping parameters

  winmin::parseandSave(argc, argv, cmd, parameters);   

  //Build the sketch for reference
  winmin::WinMin obj(parameters);
}
