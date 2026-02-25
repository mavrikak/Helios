/**
 * @file globals.h
 * @brief Global numerical constants, fundamental types, and simple macros.
 *
 * This header is included broadly across the project to provide:
 * - Fundamental **typedefs** for complex numbers and 3‑vectors.
 * - Physical **constants** (pi, i, c0, epsilon0, mu0, Z0).
 * - Minimal geometric **tags** (FRONT/BACK) used for triangle orientation.
 *
 * Keep this header lean: do **not** put run‑time configuration, wavelengths,
 * or any domain‑specific variables here.
 */

// Only define these once during compiling
#ifndef GLOBALS_H
#define GLOBALS_H

#include <blitz/array.h>
#include <blitz/tinyvec2.h>
#include <blitz/tvecglobs.h>
#include <complex>

// -----------------------------
// Types
// -----------------------------
/// Complex scalar type (double precision)
typedef std::complex<double> dcmplx;
/// Unsigned integer shorthand
typedef unsigned int uint;
/// 3‑vector of real doubles
typedef blitz::TinyVector<double, 3> rvec;
/// 3‑vector of complex doubles
typedef blitz::TinyVector<dcmplx, 3> cvec;
/// 3x3 dyadic composed of complex 3‑vectors
typedef blitz::TinyVector<cvec, 3> cdyad;

// -----------------------------
// Physical constants (SI)
// -----------------------------
/// pi = 3.141592653589793...
extern const double PI;
/// Imaginary unit @f$ i=\sqrt{-1} @f$.
extern const dcmplx I;
/// Speed of light in vacuum
extern const double CVAC;
/// Vacuum permittivity
extern const double EPS0;
/// Vacuum permeability
extern const double MU0;
/// Free‑space impedance
extern const double Z0;

// -----------------------------
// Geometric tags
// -----------------------------
/// Orientation tag for the front face of a triangle
#define FRONT 0
/// Orientation tag for the back face of a triangle
#define BACK 1

#endif
