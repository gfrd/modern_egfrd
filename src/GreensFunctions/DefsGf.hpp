#pragma once

// --------------------------------------------------------------------------------------------------------------------------------
// common includes (limit to those that are really handy)

#include <vector>
#include <string>
#include <assert.h>
#include "gfrd_compat.hpp"

// --------------------------------------------------------------------------------------------------------------------------------
// common types (prefer c++ using expression above old style typedefs)

using uint = unsigned int;
using DoubleVector = std::vector<double>;
using DoublePair = std::pair<double, double>;

// --------------------------------------------------------------------------------------------------------------------------------
// common macros (use macros definitions sparsely, its old style )

#ifndef THROW_UNLESS
#define THROW_UNLESS( CLASS, EXPRESSION ) if (!(EXPRESSION)) throw CLASS ("Check [" + std::string(#EXPRESSION) + "] failed.")
#endif

#ifndef THROW_UNLESS_MSG
#define THROW_UNLESS_MSG( CLASS, EXPRESSION , MESSAGE) if (!(EXPRESSION)) throw CLASS ("Check [" + std::string(#EXPRESSION) + "] failed: " + std::string(#MESSAGE))
#endif

#ifndef THROW_EXCEPTION
#define THROW_EXCEPTION( CLASS, MESSAGE ) { throw CLASS ("Failed: " + std::string(#MESSAGE)); }
#endif

#ifndef ASSERT
#define ASSERT( EXPRESSION )     assert( EXPRESSION );
#endif

#ifndef VERIFY
#define VERIFY( EXPRESSION )     { bool check = EXPRESSION; assert( check ); }
#endif

#ifndef UNUSED
#define UNUSED( EXPRESSION )      do { (void)(EXPRESSION); } while (0)
#endif

// --------------------------------------------------------------------------------------------------------------------------------
// common constants 

static struct GlobalGfConfig
{
   constexpr uint MAX_ORDER() const noexcept { return 50; }                      // Terms and Order are used inconsistently, need revision
   constexpr uint MIN_TERMS() const noexcept { return 20; }                      // 3DRadInf used to have MAX_ORDER=70
   constexpr uint MAX_TERMS() const noexcept { return 500; }

   const double EPSILON = 1E-12;                 // Relative numeric error (theoretical for double(64b): 2e?53 = 1.11e-16)
   const double TOLERANCE = 1E-8;                 // Error tolerance used by default.

   const double THETA_TOLERANCE = 1E-5;         // SphericalBesselGenerator's accuracy, used by some theta-related calculations
   const double MIN_T_FACTOR = 1E-8;

   const double L_TYPICAL_1D = 1E-8;             // Typical length scale of the system, may not be true!
   const double L_TYPICAL_2D = 1E-7;
   const double T_TYPICAL_1D = 1E-6;             // Typical time scale
   const double T_TYPICAL_2D = 1E-5;

   const double PDENS_TYPICAL = 1;               // Is 1E3 a good measure for the probability density?!
   const double CUTOFF_H = 6.0;                  // Cutoff distance: When H * sqrt(2Dt) < 1/2*L, use free greens function instead of absorbing.
   const double CUTOFF = 1E-10;                  // H = 4.0: ~3e-5, 4.26: ~1e-6, 5.0: ~3e-7, 5.2: ~1e-7, 5.6: ~1e-8, 6.0: ~1e-9

   const bool NO_BESSEL_TABLE = false;           // don't use precalculated bessel functions (spherical and cylindrical) but compute values with GSL (much slower)

} GfCfg;

// --------------------------------------------------------------------------------------------------------------------------------
// Library export/import definition


#if defined(_MSC_VER)
#if defined(GreensFunctions_EXPORTS)
#define GF_EXPORT __declspec(dllexport)
#else
#define GF_EXPORT __declspec(dllimport)
#endif
#else
#define GF_EXPORT
#endif

// --------------------------------------------------------------------------------------------------------------------------------
