#ifndef DEFSEGRFD_HPP
#define DEFSEGRFD_HPP

// --------------------------------------------------------------------------------------------------------------------------------
// common includes (limit to those that are really handy)

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
   const double SINGLE_SHELL_FACTOR = 3.5;
   const double MULTI_SHELL_FACTOR = std::sqrt(3);
   const double TIME_TOLERANCE = 1e-10;

   const double DEFAULT_STEPSIZE_FACTOR = 0.05;        // The maximum step size in the newBD algorithm is determined as DSSF * sigma_min.
                                                   // Make sure that DEFAULT_STEP_SIZE_FACTOR < MULTI_SHELL_FACTOR, or else the
                                                   // reaction volume sticks out of the multi.
   const double BD_DT_HARDCORE_MIN = 1e-9;         // This is to define a hardcore lower bound for the timestep that will be
                                                   // dynamically determined by the new BD scheme.It will prevent the algorithm
                                                   // to calculate ridiculously small timesteps, but will break detail balance.
                                                   // Take care : This is for testing only!Keep this at a negative value for normal sims!


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


#endif // DEFSEGRFD_HPP
