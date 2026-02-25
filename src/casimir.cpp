/**
 * @file casimir.cpp
 * @brief Main entry point for the HELIOS Casimir force routine.
 *
 * This executable drives the computation of Casimir forces via the
 * surface integral equation (SIE) solver infrastructure. Command line
 * arguments control the input job file, output prefix, solver variant,
 * accuracy settings, and threading configuration.
 */

#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include "CLParser.h"
#include "GreenF.h"
#include "JobParser.h"
#include "SIEFormPMCHW.h"
#include "SimJob.h"
#include "SurfaceMesh.h"

/**
 * @brief Main driver routine for the simulation program.
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return 0 on success.
 *
 * Parses command line arguments, determines job type, and runs the requested
 * simulation or mesh conversion. Provides help text when requested.
 *
 * Supported command line switches:
 * -j FILE : Specify job file (default: input.job)
 * -op PREFIX: Output prefix for result files
 * -l [0,1,2]: Solver type (0=direct, 1=cgsqr, 2=LU [default])
 * -t FILE : Look-up table for error functions (optional)
 * -a : Enable high-accuracy integration for matrix/vector elements
 * -etm N : Number of evaluation terms in periodic simulation (default 1)
 * -th N : Number of threads (0=all cores, default=1)
 * --help : Print usage message
 */
int main(int argc, char *argv[]) {
  // ----------------------
  // Defaults for arguments
  // ----------------------
  std::string jobFile("input.job"); //!< Input job specification
  std::string outPrefix("");        //!< Prefix for output files
  bool help(false);                 //!< Help flag
  bool NeedsAccurate(false);        //!< Request high-accuracy quadrature
  int solver(2);                    //!< Solver index (default: LU=2)
  int etm(1);                       //!< Number of periodic evaluation terms
  std::string tableFile("");        //!< Lookup table filename (optional)
  int threads(1);                   //!< Number of parallel threads

  // ----------------------
  // Declare supported CLI options
  // ----------------------
  std::vector<CLParser *> clParser;
  clParser.push_back(new CLStringOption("-j", &jobFile));
  clParser.push_back(new CLStringOption("-op", &outPrefix));
  clParser.push_back(new CLOption<int>("-l", &solver));
  clParser.push_back(new CLStringOption("-t", &tableFile));
  clParser.push_back(new CLFlag("--help", &help));
  clParser.push_back(new CLFlag("-a", &NeedsAccurate));
  clParser.push_back(new CLOption<int>("-etm", &etm));
  clParser.push_back(new CLOption<int>("-th", &threads));

  // ----------------------
  // Parse CLI and override defaults
  // ----------------------
  for (std::vector<CLParser *>::iterator pIter = clParser.begin();
       pIter != clParser.end(); pIter++) {
    (*pIter)->Parse(argc, argv);
  }

  // ----------------------
  // Help or execution
  // ----------------------
  if (help) {
    std::cout << "HELIOS SIE program : Casimir computation" << std::endl;
  } else {
    SimJob *job = new SimJob;

    // Enable more accurate quadrature if requested
    if (NeedsAccurate) GreenF::EnableAccurate();

    // Sanitize thread count
    if (threads < 1 or threads > std::thread::hardware_concurrency()) {
      threads = std::thread::hardware_concurrency();
    }
    SIEFormPMCHW::AssignThreads(threads);
    SimJob::AssignThreads(threads);

    // Run Casimir computation
    return job->Casimir(jobFile, outPrefix);
  }
  return 0;
}

