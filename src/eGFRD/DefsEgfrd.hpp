#pragma once

// --------------------------------------------------------------------------------------------------------------------------------
// common includes (limit to those that are really handy)

#include <iostream>
#include <cassert>
#include <cmath>
#include <csignal>
#include "exceptions.hpp"
#include "makeString.hpp"
#include "gfrd_compat.hpp"

// --------------------------------------------------------------------------------------------------------------------------------
// common types (prefer c++ using expression above old style typedefs)

using uint = unsigned int;
using idtype = unsigned long long;

// --------------------------------------------------------------------------------------------------------------------------------
// common macros (use macros definitions sparsely, its old style )

#ifndef THROW_UNLESS
#define THROW_UNLESS( CLASS, EXPRESSION ) if (!(EXPRESSION)) { DEBUGBREAK; std::string ex_msg = make_string() << "EXCEPTION: " << #CLASS << " in " << __FILE__ << ":" << __LINE__ <<  " - " << __FUNCTION__ << "() " << "Check [" << #EXPRESSION << "] failed."; throw CLASS (ex_msg); } 
#endif

#ifndef THROW_UNLESS_MSG
#define THROW_UNLESS_MSG( CLASS, EXPRESSION , MESSAGE) if (!(EXPRESSION)) { DEBUGBREAK; std::string ex_msg = make_string() << "EXCEPTION: " << #CLASS << " in " << __FILE__ << ":" << __LINE__ <<  " - " << __FUNCTION__ << "() " << MESSAGE; throw CLASS (ex_msg); }
#endif

#ifndef THROW_EXCEPTION
#define THROW_EXCEPTION( CLASS, MESSAGE ) { DEBUGBREAK; std::string ex_msg = make_string() << "EXCEPTION: " << #CLASS  << " in " << __FILE__ << ":" << __LINE__ <<  " - " << __FUNCTION__ << "() " << MESSAGE; throw CLASS (ex_msg); }
#endif

#ifndef ASSERT
#define ASSERT( EXPRESSION )     assert( EXPRESSION );
#endif

#ifndef VERIFY
#define VERIFY( EXPRESSION )     { bool check = EXPRESSION; assert( check ); UNUSED(check) }
#endif

#ifndef UNUSED
#define UNUSED( ... )     { unused(__VA_ARGS__); }
#endif

#ifndef DEBUGBREAK
#ifdef _DEBUG
#ifdef _MSC_VER
#define DEBUGBREAK         __debugbreak();
#else
#define DEBUGBREAK         std::raise(SIGTRAP);
#endif
#else
#define DEBUGBREAK
#endif
#endif

#ifndef M_PI
#define M_PI       3.14159265358979323846264338328      /* pi */
#endif


template<class... T> void unused(T&&...) { }

// --------------------------------------------------------------------------------------------------------------------------------
// common constants 

static struct GlobalGfrdConfig
{
   const double MINIMAL_SEPARATION_FACTOR = 1.0 + 1e-07;
   const double SAFETY = 1.0 + 1e-2;
   const double SINGLE_SHELL_FACTOR = 3.5;         // This is the threshold for when the algorithm switches from forming
                                                   // NonInteractionSingles to forming a Pair or Interaction. It also defines
                                                   // the radius in which the NonInteractionSingle will burst intruding domains.
                                                   // IMPORTANT NOTE: SINGLE_SHELL_FACTOR should be AT_LEAST 2 * MULTI_SHELL_FACTOR !
                                                   // Otherwise we risk the construction of situations in which two neighbouring
                                                   // particles neither can burst nor form a multi nor form a minimal single!

   const double MULTI_SHELL_FACTOR = std::sqrt(3); // This factor multiplied with the particle radius decides when to add
                                                   // NonInteractionSingles to a Multi and also defines the Multi shell size.
                                                   // IMPORTANT NOTE: MULTI_SHELL_FACTOR should be AT LEAST sqrt(2) !
                                                   // This stems from the fact that there is vacant space in the cylinder //
   const double TIME_TOLERANCE = 1e-10;
   const double MAX_CELL_OCCUPANCY = 0.5;          // Maximum size of domains relative to world cell size

   const double DEFAULT_STEPSIZE_FACTOR = 0.05;    // The maximum step size in the newBD algorithm is determined as DSSF * sigma_min.
                                                   // Make sure that DEFAULT_STEP_SIZE_FACTOR < MULTI_SHELL_FACTOR, or else the
                                                   // reaction volume sticks out of the multi.
   const double BD_DT_HARDCORE_MIN = -1e-9;        // This is to define a hardcore lower bound for the timestep that will be
                                                   // dynamically determined by the new BD scheme. It will prevent the algorithm
                                                   // to calculate ridiculously small timesteps, but will break detail balance.
                                                   // Take care : This is for testing only! Keep this at a negative value for normal sims!
   const uint ThrowInRetryCount = 500000;          // When randomly inserting particles, stop after N failed placements per particle (due to overlap or out-of-bounds error)

} GfrdCfg;

// --------------------------------------------------------------------------------------------------------------------------------
// Library export/import definition

#if defined(_MSC_VER)
#if defined(eGFRD_EXPORTS)
#define GFRD_EXPORT __declspec(dllexport)
#else
#define GFRD_EXPORT __declspec(dllimport)
#endif
#else
#define GFRD_EXPORT
#endif

// --------------------------------------------------------------------------------------------------------------------------------

// forward declare available boundary conditions
struct WorldCyclic;
struct WorldNoBounds;

// Specify customized Matrix Cell grid on build command line (e.q. make -DMATRIXSIZE=4)
#ifndef MATRIXSIZE
#define MATRIXSIZE 8
#endif

// Specify build configuration of world/simulator types (changing these needs rebuild of sources, only one kind of simulator per build is supported)
struct CompileConfigSimulator
{
   using TBoundCondition = WorldCyclic;
   static const size_t MatrixCellsX = MATRIXSIZE;
   static const size_t MatrixCellsY = MATRIXSIZE;
   static const size_t MatrixCellsZ = MATRIXSIZE;
};

// --------------------------------------------------------------------------------------------------------------------------------

static_assert(CompileConfigSimulator::MatrixCellsX >= 3, "MatrixSpace cells should be equal or larger than 3");      // in a cyclic world smaller matrixspace is no optimization
static_assert(CompileConfigSimulator::MatrixCellsY >= 3, "MatrixSpace cells should be equal or larger than 3");      // because you need to check all 26 neighbor spaces anyway
static_assert(CompileConfigSimulator::MatrixCellsZ >= 3, "MatrixSpace cells should be equal or larger than 3");      // use power of 2 for best performance

// --------------------------------------------------------------------------------------------------------------------------------

class Persistence;

// --------------------------------------------------------------------------------------------------------------------------------

inline void gfrd_print_header()
{
   std::cout << "eGFRD (modern), NWO-I AMOLF, 2017, Amsterdam, The Netherlands, www.grfd.org" << std::endl;
}

inline void gfrd_print_version()
{
   gfrd_print_header();
   std::cout << " version: " << GFRD_VERSION_MAJOR << "." << GFRD_VERSION_MINOR << "." << GFRD_VERSION_PATCH << std::endl;
   std::cout << " build: " << GFRD_VERSION_BUILD << std::endl;
   std::cout << " platform: " << (sizeof(std::nullptr_t) == 8 ? "x64" : "x86");
#if _DEBUG
   std::cout << " DEBUG";
#endif
   std::cout << std::endl << " matrixspace: " << CompileConfigSimulator::MatrixCellsX << "x" << CompileConfigSimulator::MatrixCellsY << "x" << CompileConfigSimulator::MatrixCellsZ << std::endl;
   std::cout << std::endl;
}

// --------------------------------------------------------------------------------------------------------------------------------
