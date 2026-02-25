/**
 * @file maintable.cpp
 * @brief Command-line utility to generate, read, and write periodic error-function
 * lookup tables used in the SIE solver. Supports both text and binary
 * formats for tables. Uses mathFunctions.h (LookupTable classes and
 * ErfcLookupTable generators).
 */

#include <blitz/array.h>
#include <iostream>
#include <string>
#include <cmath>
#include <stdio.h>
#include "CLParser.h"
#include "mathFunctions.h"

/**
 * @brief Main driver for periodic error-function lookup table management.
 * @param argc Argument count.
 * @param argv Argument vector.
 * @return 0 on success.
 *
 * Supports the following modes (mutually exclusive):
 * --help : Print usage message
 * -r FILE : Read an existing text lookup table
 * -w FILE : Write a new text lookup table (user prompted for bounds)
 * -rb FILE : Read an existing binary lookup table
 * -wb FILE : Write a new binary lookup table (user prompted for bounds)
 *
 * Text and binary formats store sampled values of the error function over
 * a 2D grid in the complex plane (Re, Im).
 */
int main(int argc, char *argv[]) {
	bool help(false);
	std::string fileName1("");
	std::string fileName2("");
	std::string fileName3("");
	std::string fileName4("");

	// -------------------------------------------------------------------------
	// Setup command line arguments
	// -------------------------------------------------------------------------
	std::vector<CLParser *> clParser;
	clParser.push_back(new CLStringOption("-r", &fileName1));
	clParser.push_back(new CLStringOption("-w", &fileName2));
	clParser.push_back(new CLStringOption("-rb", &fileName3));
	clParser.push_back(new CLStringOption("-wb", &fileName4));
	clParser.push_back(new CLFlag("--help", &help));

	// Parse command line and overwrite defaults
	for (std::vector<CLParser *>::iterator pIter = clParser.begin();
			pIter != clParser.end(); pIter++) {
		(*pIter)->Parse(argc, argv);
	}

	// -------------------------------------------------------------------------
	// Dispatch modes
	// -------------------------------------------------------------------------
	if (help) {
		std::cout
				<< "HELIOS SIE periodic error function table generation"
				<< std::endl;
		std::cout << "   read text table: -r [FILENAME]" << std::endl;
		std::cout << "   write text table: -w [FILENAME]" << std::endl;
		std::cout << "   read binary table: -rb [FILENAME]" << std::endl;
		std::cout << "   write binary table: -wb [FILENAME]" << std::endl;
	} else if (fileName1 != "") {
		std::cout << "Loading lookup table..." << std::endl;
		LookupTable* t = new LookupTable(fileName1);
		std::cout << "Size:" << t->sizeRe << "," << t->sizeIm << std::endl;
		std::cout << "Max:" << t->maxRe << "," << t->maxIm << std::endl;
		std::cout << "Inc:" << t->incRe << "," << t->incIm << std::endl;
		std::cout << "f" << t->arg[0][0] << ":" << t->fun[0][0] << std::endl;

	} else if (fileName2 != "") {
		double maxRe, maxIm, incRe, incIm;
		std::cout << "Enter maxRe maxIm incRe incIm:" << std::endl;
		std::cin >> maxRe >> maxIm >> incRe >> incIm;

		ErfcLookupTable(fileName2, maxRe, maxIm, incRe, incIm);

	} else if (fileName3 != "") {
		std::cout << "Loading lookup table..." << std::endl;
		LookupTableBin* tb = new LookupTableBin(fileName3);
		std::cout << "Size:" << tb->sizeRe << "," << tb->sizeIm << std::endl;
		std::cout << "Max:" << tb->maxRe << "," << tb->maxIm << std::endl;
		std::cout << "Inc:" << tb->incRe << "," << tb->incIm << std::endl;
		std::cout << "f" << tb->arg[0][0] << ":" << tb->fun[0][0] << std::endl;

	} else if (fileName4 != "") {
		double maxRe, maxIm, incRe, incIm;
		std::cout << "Enter maxRe maxIm incRe incIm:" << std::endl;
		std::cin >> maxRe >> maxIm >> incRe >> incIm;
		ErfcLookupTableBin(fileName4, maxRe, maxIm, incRe, incIm);
	}

}

