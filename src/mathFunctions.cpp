#include "mathFunctions.h"
#include <blitz/tvecglobs.h>
#include <fstream>
#include <iostream>
#include "iofunctions.h"

// -----------------------------------------------------------------------------
// SymmetryOperation implementations
// -----------------------------------------------------------------------------
Translation::Translation(rvec translationVector) : t(translationVector){};

rvec Translation::Apply(rvec inputVector) {
  rvec out(inputVector);
  out -= t;
  return out;
}

// -----------------------------------------------------------------------------
// Special function implementations
// -----------------------------------------------------------------------------
// Error function erf(z) for complex arguments (series and asymptotics).
dcmplx cerf(dcmplx arg) {
  dcmplx c0, z1, result;
  c0 = exp(-arg * arg);

  if (arg.real() < 0)
    z1 = -arg;
  else
    z1 = arg;

  if (abs(arg) < 5.8) {
    dcmplx cs = z1;
    dcmplx cr = z1;
    for (int k = 0; k < 120; k++) {
      cr *= z1 * z1 / (k + 1.5);
      cs += cr;
      if (abs(cr / cs) < 1e-15) break;
    }
    result = 2. * c0 * cs / sqrt(PI);
  } else {
    dcmplx cl = 1. / z1;
    dcmplx cr = cl;
    for (int k = 0; k < 13; k++) {
      cr *= -(k + 0.5) / (z1 * z1);
      cl += cr;
      if (abs(cr / cl) < 1e-15) break;
    }
    result = 1. - c0 * cl / sqrt(PI);
  }

  if (arg.real() < 0) return -result;

  return result;
}

// Complementary error function erfc(z).
dcmplx cerfc(dcmplx arg) {
  dcmplx id(1., 0.);
  return id - cerf(arg);
}

// Exponential integral Ei(x) for positive real x.
double Ei(const double &x) {
  double EI_TOLERANCE = 1e-9;
  double EULER = 0.57721566490153;
  int n;
  double Ei, term;
  Ei = 0.0;
  if (x < 15.0) {
    n = 1;
    term = x;
    Ei = EULER + log(x) + term;
    while (fabs(term) > EI_TOLERANCE) {
      term *= x * n / (n + 1) / (n + 1);
      Ei += term;
      n++;
    }
  }
  return Ei;
}

// Exponential integral Ep(x,p) of order p.
dcmplx Ep(const double &x, int p) {
  int n(1);
  dcmplx Ep;
  if (x >= 0) {
    double E1_TOLERANCE = 1e-9;
    double EULER = 0.57721566490153;
    int m;
    double E1, term;
    E1 = 0.0;
    if (x < 15.0) {
      m = 1;
      term = -x;
      E1 = -EULER - log(x) - term;
      while (fabs(term) > E1_TOLERANCE) {
        term *= -x * m / (m + 1) / (m + 1);
        E1 -= term;
        m++;
      }
    }
    Ep = E1;
  } else {
    Ep = -Ei(-x) + I * PI;
  }
  while (n < p) {
    Ep = 1. / n * (exp(-x) - x * Ep);
    ++n;
  }
  return Ep;
}

// Compute factorial of integer n.
int factorial(int num) {
  int result(1);
  for (int i = 1; i <= num; ++i) result *= i;
  return result;
}

// Transpose of a complex dyadic (3x3 matrix).
cdyad transpose(cdyad d) {
  return cdyad(cvec(d(0)(0), d(1)(0), d(2)(0)), cvec(d(0)(1), d(1)(1), d(2)(1)),
               cvec(d(0)(2), d(1)(2), d(2)(2)));
}

// -----------------------------------------------------------------------------
// Table generation routines
// -----------------------------------------------------------------------------
// Generate and save text lookup table of cerfc values.
int ErfcLookupTable(std::string filename, double maxRe, double maxIm,
                    double incRe, double incIm) {
  dcmplx z;
  std::vector<std::vector<dcmplx> > erfcTable;
  std::vector<std::vector<dcmplx> > posTable;
  std::cout << "Computing look-up table" << std::endl;
  for (double x = -maxRe; x <= maxRe; x += incRe) {
    std::cout << "." << std::flush;
    std::vector<dcmplx> buffer, buffer1;
    for (double y = -maxIm; y <= maxIm; y += incIm) {
      z = x + I * y;
      buffer.push_back(cerfc(z));
      buffer1.push_back(z);
    }
    erfcTable.push_back(buffer);
    posTable.push_back(buffer1);
  }
  std::ofstream outputFile(filename.c_str());
  std::cout << "finished...writing file " << filename << std::endl;
  if (!outputFile.is_open()) {
    std::cout << "Error opening erfc output file: " << filename << std::endl;
    return 1;
  }
  outputFile << "# Look-up table for erfc function" << std::endl;
  outputFile << "# Table parameters: size along x, size along y" << std::endl;
  outputFile << erfcTable.size() << " " << erfcTable[0].size() << std::endl;
  outputFile << "# Table: x, y, re(erfc(x+iy)),im(erfc(x+iy))" << std::endl;
  for (unsigned int i(0); i < erfcTable.size(); ++i) {
    for (unsigned int j(0); j < erfcTable[i].size(); ++j) {
      outputFile << posTable[i][j] << " " << erfcTable[i][j] << std::endl;
    }
  }
  outputFile.close();
  return 0;
}

// Construct LookupTable from text file.
LookupTable::LookupTable()
    : sizeRe(0), sizeIm(0), maxRe(0), maxIm(0), incRe(0), incIm(0) {}

// Construct LookupTable from text file.
LookupTable::LookupTable(std::string fileName) {
  std::ifstream fin(fileName.c_str());
  if (!fin.is_open()) {
    std::cout << "Error opening look-up table file: " << fileName << std::endl;
  }
  dcmplx z;
  int err = ReadCommented<int>(&fin, &sizeRe);
  err += ReadCommented<int>(&fin, &sizeIm);
  for (int i(0); i < sizeRe; ++i) {
    std::vector<dcmplx> buffer, buffer1;
    for (int j(0); j < sizeIm; ++j) {
      err += ReadCommented<dcmplx>(&fin, &z);
      buffer1.push_back(z);
      err += ReadCommented<dcmplx>(&fin, &z);
      buffer.push_back(z);
    }
    fun.push_back(buffer);
    arg.push_back(buffer1);
  }
  maxRe = -real(arg[0][0]);
  maxIm = -imag(arg[0][0]);
  incRe = real(arg[1][0] - arg[0][0]);
  incIm = imag(arg[0][1] - arg[0][0]);
  fin.close();
}

// Generate and save binary lookup table of cerfc values.
int ErfcLookupTableBin(std::string filename, double maxRe, double maxIm,
                    double incRe, double incIm) {
  dcmplx z;
  std::vector<std::vector<dcmplx> > erfcTable;
  std::vector<std::vector<dcmplx> > posTable;
  std::cout << "Computing look-up table" << std::endl;
  for (double x = -maxRe; x <= maxRe; x += incRe) {
    std::cout << "." << std::flush;
    std::vector<dcmplx> buffer, buffer1;
    for (double y = -maxIm; y <= maxIm; y += incIm) {
      z = x + I * y;
      buffer.push_back(cerfc(z));
      buffer1.push_back(z);
    }
    erfcTable.push_back(buffer);
    posTable.push_back(buffer1);
  }
  std::cout << "finished...writing file " << filename << std::endl;

  FILE* fp = fopen(filename.c_str(), "wb");
  if (fp==NULL) {
    std::cout << "Error opening erfc output file: " << filename << std::endl;
    return 1;
  }
  int temp = erfcTable.size();
  fwrite(&temp, sizeof(int), 1, fp);
  temp = erfcTable[0].size();
  fwrite(&temp, sizeof(int), 1, fp);
  double temp2;
  for (unsigned int i(0); i < erfcTable.size(); ++i) {
    for (unsigned int j(0); j < erfcTable[i].size(); ++j) {
      temp2 = posTable[i][j].real();
      fwrite(&temp2, sizeof(double), 1, fp);
      temp2 = posTable[i][j].imag();
      fwrite(&temp2, sizeof(double), 1, fp);
      temp2 = erfcTable[i][j].real();
      fwrite(&temp2, sizeof(double), 1, fp);
      temp2 = erfcTable[i][j].imag();
      fwrite(&temp2, sizeof(double), 1, fp);
    }
  }
  fclose(fp);
  return 0;
}

// Construct LookupTable from bin file.
LookupTableBin::LookupTableBin()
    : sizeRe(0), sizeIm(0), maxRe(0), maxIm(0), incRe(0), incIm(0) {}

// Construct LookupTable from bin file.
LookupTableBin::LookupTableBin(std::string fileName) {

  FILE* fp = fopen(fileName.c_str(), "rb");

  double temp1,temp2;
  dcmplx z;

  if (fp==NULL) {
    std::cout << "Error opening look-up table file: " << fileName << std::endl;
  }

  int err = fread(&sizeRe, sizeof(int), 1, fp);
  err += fread(&sizeIm, sizeof(int), 1, fp);

  for (int i(0); i < sizeRe; ++i) {
    std::vector<dcmplx> buffer, buffer1;
    for (int j(0); j < sizeIm; ++j) {
      err += fread(&temp1, sizeof(double), 1, fp);
      err += fread(&temp2, sizeof(double), 1, fp);
      z = temp1 + I * temp2;
      buffer1.push_back(z);
      err += fread(&temp1, sizeof(double), 1, fp);
      err += fread(&temp2, sizeof(double), 1, fp);
      z = temp1 + I * temp2;
      buffer.push_back(z);
    }
    fun.push_back(buffer);
    arg.push_back(buffer1);
  }
  maxRe = -real(arg[0][0]);
  maxIm = -imag(arg[0][0]);
  incRe = real(arg[1][0] - arg[0][0]);
  incIm = imag(arg[0][1] - arg[0][0]);
  fclose(fp);
}
