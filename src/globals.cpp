/**
 * @file globals.cpp
 * @brief Definitions for global constants declared in globals.h.
 */

#include "globals.h"

// Numerical constants (SI)
/** 
* @brief Imaginary unit. 
* @return Imaginary unit @f$ i=\sqrt{-1} @f$.
*/
const dcmplx I(0., 1.);                                     ///< Imaginary unit i
const double PI = 3.14159265358979323846264338328;          ///< pi
const double CVAC = 299792458.;                             ///< Speed of light in vacuum
const double EPS0 = 8.85418781762038985053656303171e-12;    ///< Vacuum permittivity
const double MU0 = 1.25663706143591729538505735331e-6;      ///< Vacuum permeability
const double Z0 = 376.730313461770655468198400420;          ///< Free-space impedance
