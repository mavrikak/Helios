/**
 * @file mainpp.cpp
 * @brief This is the main file for the HELIOS program,
 * specifically the post-processing portion. Handles field evaluation 
 * after a solution has been computed. Allows specification of which 
 * quantities to output.
 */

#include <iostream>
#include <string>
#include <thread>
#include <vector>
#include "CLParser.h"
#include "SIEFormPMCHW.h"
#include "SimJob.h"

/**
 * @brief Main driver for post-processing in the HELIOS SIE program.
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return 0 on success.
 *
 * This executable evaluates electromagnetic fields after a simulation job
 * has been completed, using the stored solution and specified observation
 * points. Users can select which field components to compute and output. Also, 
 * it evaluates other quantities such as power and cross sections.
 *
 * Command line switches:
 * -j FILE : Job file (default input.job)
 * -s FILE : Solution file (default output.sol)
 * -p FILE : Positions file (default points.pos)
 * -op PREFIX: Output prefix for result files
 * -t FILE : Look-up table file (optional)
 * -a : Enable high-accuracy integration
 * -etm N : Evaluation terms for periodic sim (default 1)
 * -th N : Number of threads (default 1, 0=all cores)
 * -fE/-fH : Compute electric/magnetic fields
 * -fD/-fB : Compute displacement/flux density fields
 * -fP : Compute Poynting vector
 * -fPw : Compute power crossing domain boundaries
 * --help : Print usage
 */
int main(int argc, char *argv[]) {
  // -------------------------------------------------------------------------
  // Defaults for job configuration
  // -------------------------------------------------------------------------
  std::string jobFile("input.job");
  std::string solutionFile("output.sol");
  std::string positionFile("points.pos");
  std::string outPrefix("");
  std::string tableFile("");
  bool help(false);
  bool NeedsAccurate(false);
  int etm(1);
  int threads(1);

  // Flags selecting which field types to compute
  bool fldE(false);
  bool fldH(false);
  bool fldD(false);
  bool fldB(false);
  bool fldP(false);
  bool fldPw(false);

  // -------------------------------------------------------------------------
  // Setup command line argument parsers
  // -------------------------------------------------------------------------
  std::vector<CLParser *> clParser;
  clParser.push_back(new CLStringOption("-j", &jobFile));
  clParser.push_back(new CLStringOption("-s", &solutionFile));
  clParser.push_back(new CLStringOption("-p", &positionFile));
  clParser.push_back(new CLStringOption("-op", &outPrefix));
  clParser.push_back(new CLFlag("--help", &help));
  clParser.push_back(new CLFlag("-fE", &fldE));
  clParser.push_back(new CLFlag("-fH", &fldH));
  clParser.push_back(new CLFlag("-fD", &fldD));
  clParser.push_back(new CLFlag("-fB", &fldB));
  clParser.push_back(new CLFlag("-fP", &fldP));
  clParser.push_back(new CLFlag("-fPw", &fldPw));
  clParser.push_back(new CLStringOption("-t", &tableFile));
  clParser.push_back(new CLFlag("-a", &NeedsAccurate));
  clParser.push_back(new CLOption<int>("-etm", &etm));
  clParser.push_back(new CLOption<int>("-th", &threads));

  // Parse command line
  for (std::vector<CLParser *>::iterator pIter = clParser.begin();
       pIter != clParser.end(); pIter++) {
    (*pIter)->Parse(argc, argv);
  }

  // -------------------------------------------------------------------------
  // Build field mask bitset
  // -------------------------------------------------------------------------
  uint fields(0);
  if ((fldE || fldH || fldD || fldB || fldP || fldPw) == false)
    // Defaults: E and H fields
    fields = EFIELD | HFIELD;
  else {
    // If field flags are given:
    if (fldE) fields |= EFIELD;
    if (fldH) fields |= HFIELD;
    if (fldD) fields |= DFIELD;
    if (fldB) fields |= BFIELD;
    if (fldP) fields |= POYVEC;
    if (fldPw) fields |= POWER;
  }

  std::cout << "Output fields type: " << fields << std::endl;
  
  // -------------------------------------------------------------------------
  // Help or execution
  // -------------------------------------------------------------------------
  if (help) {
    std::cout << "HELIOS SIE simulation program, field output"
              << std::endl;
    std::cout << "   usage: mainpp [-j JOBFILE] [-p POSITIONFILE] [fields]"
              << std::endl;
    std::cout << "   defaults:      -j input.job -p points.pos    -fE -fH"
              << std::endl;
    std::cout << "   possible fields are E, H, D, B, P (Poynting vectors)"
              << std::endl;
    std::cout << "   -fPw : power crossing the boundary of the domains"
              << std::endl;
    std::cout << "   switches:    -op PREFIX" << std::endl;
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
    std::cout << "                -th [1,2,3,...]" << std::endl;
    std::cout << "                   number of parallel threads to run, "
                 "default value is 1"
              << std::endl;
  } else {
    // -----------------------------------------------------------------------
    // Launch field evaluation
    // -----------------------------------------------------------------------
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
    return job->FieldEval(jobFile, positionFile, outPrefix, fields);
  }
  return 0;
}

