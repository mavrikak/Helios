/**
 * @file main.cpp
 * @brief This is the main file for the HELIOS program. 
 * Provides the entry point to the SIE (Surface Integral Equation) solver.
 * Handles command-line argument parsing, dispatches jobs, and manages
 * post-processing routines.
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
 * @brief Main driver routine for the HELIOS SIE simulation program.
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
  // -------------------------------------------------------------------------
  // Defaults for job configuration
  // -------------------------------------------------------------------------
  std::string jobFile("input.job"); ///< Input job filename
  std::string outPrefix("");        ///< Prefix for output files
  bool help(false);                 ///< Help flag
  bool NeedsAccurate(false);        ///< Flag for high-accuracy evaluation
  int solver(2);                    ///< Solver type (default LU=2)
  int etm(1);                       ///< Number of evaluation terms (periodic)
  std::string tableFile("");        ///< Optional look-up table file
  int threads(1);                   ///< Number of threads

  // -------------------------------------------------------------------------
  // Setup command line parsers
  // -------------------------------------------------------------------------
  std::vector<CLParser *> clParser;
  clParser.push_back(new CLStringOption("-j", &jobFile));
  clParser.push_back(new CLStringOption("-op", &outPrefix));
  clParser.push_back(new CLOption<int>("-l", &solver));
  clParser.push_back(new CLStringOption("-t", &tableFile));
  clParser.push_back(new CLFlag("--help", &help));
  clParser.push_back(new CLFlag("-a", &NeedsAccurate));
  clParser.push_back(new CLOption<int>("-etm", &etm));
  clParser.push_back(new CLOption<int>("-th", &threads));

  // Parse command line and overwrite defaults
  for (std::vector<CLParser *>::iterator pIter = clParser.begin();
       pIter != clParser.end(); pIter++) {
    (*pIter)->Parse(argc, argv);
  }

  // -------------------------------------------------------------------------
  // Handle help flag
  // -------------------------------------------------------------------------
  if (help) {
    std::cout << "HELIOS SIE simulation program" << std::endl;
    std::cout << "   usage: main [-j JOBFILE]" << std::endl;
    std::cout << "   defaults:    -j input.job" << std::endl;
    std::cout << "   switches:    -l [0,1,2]" << std::endl;
    std::cout << "                   specify solver type:" << std::endl;
    std::cout << "                   0: direct solver" << std::endl;
    std::cout << "                   1: cgsqr " << std::endl;
    std::cout << "                   2: LU decomposition [default]"
              << std::endl;
    std::cout << "                -op PREFIX" << std::endl;
    std::cout << "                   specify prefix for output file(s)"
              << std::endl;
    std::cout << "                -t look-up table file" << std::endl;
    std::cout << "                   look-up table for error functions - not "
                 "used by default"
              << std::endl;
    std::cout << "                -a high-accuracy integration for matrix and "
                 "vector elements"
              << std::endl;
    std::cout << "                -etm [1,2,3,...]" << std::endl;
    std::cout << "                   number of evaluation terms in periodic "
                 "simulation, default value is 1"
              << std::endl;
    std::cout << "                -th [0,1,2,...]" << std::endl;
    std::cout << "                   number of parallel threads to run, "
                 "default value is 1, 0 uses maximum threads"
              << std::endl;
  } else {
    // -----------------------------------------------------------------------
    // Parse job file to decide task
    // -----------------------------------------------------------------------
    JobParser J;
    J.open(jobFile);
    std::string JobType(J.readTag());
    J.close();
    if (J.getTag() == "conversion") { // Mesh conversion job
      SurfaceMesh s;
      return s.ConvertMeshFile(jobFile);
    } else if (J.getTag() == "job") { // Simulation job
      SimJob *job = new SimJob;
      if (tableFile != "") job->LoadTableFromFile(tableFile);
      if (NeedsAccurate) Domain::EnableAccurate();
      if (NeedsAccurate) GreenF::EnableAccurate();
      if (NeedsAccurate) IncidentField::EnableAccurate();
      if (etm < 1) {
        etm = 1;
      }
      GreenF::AssignEtm(etm);
      if (threads < 1 or threads > std::thread::hardware_concurrency()) {
        threads = std::thread::hardware_concurrency();
      }
      SIEFormPMCHW::AssignThreads(threads);
      SimJob::AssignThreads(threads);

      // Run the simulation
      if ( job->Simulate(jobFile, solver, outPrefix) == 0) {
        std::cout << "Simulation completed successfully." << std::endl;

        // Run post-processing
        SimJob *jobPP = new SimJob;
        if (tableFile != "") job->LoadTableFromFile(tableFile);
        if (NeedsAccurate) Domain::EnableAccurate();
        if (NeedsAccurate) GreenF::EnableAccurate();
        if (NeedsAccurate) IncidentField::EnableAccurate();
        if (etm < 1) {
          etm = 1;
        }
        GreenF::AssignEtm(etm);
        if (threads < 1 or threads > std::thread::hardware_concurrency()) {
          threads = std::thread::hardware_concurrency();
        }
        SIEFormPMCHW::AssignThreads(threads);
        SimJob::AssignThreads(threads);
        return jobPP->FieldEval(jobFile, "", "", 0, true);
      } else {
        std::cout << "Error during simulation." << std::endl;
        return 1;
      }
    } else {
      std::cout << "Error ! Unknown opening tag." << std::endl;
      return 1;
    }
  }
  return 0;
}

