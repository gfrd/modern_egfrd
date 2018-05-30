#include <gsl/gsl_sum.h>
#include <gsl/gsl_errno.h>
#include <cmath>
#include "DefsGf.hpp"
#include "helperFunctionsGf.hpp"
#include "Logger.hpp"
#include "SphericalBesselGenerator.hpp"
#include "CylindricalBesselGenerator.hpp"
#include <iomanip>
#include <fstream>

#if defined(_MSC_VER)
#include <windows.h>
#else
#include <unistd.h>
#include <cstring>
#define MAX_PATH 1024
#endif

// --------------------------------------------------------------------------------------------------------------------------------

std::string GetExecutablePath()
{
   char buffer[MAX_PATH];
#if defined(_MSC_VER)
   const char sep = '\\';
   auto module = GetModuleHandleA("greensfunctions.dll");
   auto ret = GetModuleFileNameA(module, buffer, MAX_PATH);
   if (ret)
#else
   const char sep = '/';
   auto ret = ::readlink("/proc/self/exe", buffer, MAX_PATH);
   if (ret != -1)
#endif
   {
      auto pos = std::strrchr(buffer, sep);
      if (pos) *(pos + 1) = '\0';
   }
   else
      buffer[0] = '\0';

   return std::string(buffer);
}

// --------------------------------------------------------------------------------------------------------------------------------

bool file_exists(const std::string& name)
{
   std::ifstream file(name);
   return file.good();	
}

// --------------------------------------------------------------------------------------------------------------------------------

double funcSum_all(std::function<double(uint i)> f, uint max_i)
{
   double p = f(0);
   if (p == 0.0) return 0.0;

   double sum = p;
   for (uint i = 1; i < max_i; ++i)
      sum += f(i);
   return sum;
}

// --------------------------------------------------------------------------------------------------------------------------------

// Will simply calculate the sum over a certain function f, until it converges
// (i.e. the sum > tolerance*current_term for a CONVERGENCE_CHECK number of 
// terms), or a maximum number of terms is summed (usually 2000).
//
// Input:
// - f: A       function object
// - max_i:     maximum number of terms it will evaluate
// - tolerance: convergence condition, default value 1-e8 (see .hpp)
double funcSum(std::function<double(uint i)> f, uint max_i, double tolerance)
{
   uint CONVERGENCE_CHECK = 4;
   double p = f(0);
   if (p == 0.0) return 0.0;

   DoubleVector table;
   table.emplace_back(p);

   double sum = p;
   bool extrapolationNeeded = true;
   uint convergenceCounter = 0;

   uint i = 1;
   for (; i < max_i;)
   {
      p = f(i);
      table.emplace_back(p);
      sum += p;
      ++i;

      if (std::fabs(sum) * tolerance >= std::fabs(p))
         convergenceCounter++;
      else
         convergenceCounter = 0;

      if (convergenceCounter >= CONVERGENCE_CHECK)
      {
         extrapolationNeeded = false;
         break;
      }
   }

   if (extrapolationNeeded)
   {
      double error;
      gsl_sum_levin_utrunc_workspace* workspace = gsl_sum_levin_utrunc_alloc(i);
      gsl_sum_levin_utrunc_accel(&table[0], table.size(), workspace, &sum, &error);
      if (std::fabs(error) >= std::fabs(sum * tolerance * 10))
         Log("funcSum").error() << "series acceleration failed! sum: " << std::scientific << std::setprecision(16) << sum << ", error: " << std::fabs(error) << " (rel error: " << std::fabs(error / sum) << "), terms_used: " << workspace->terms_used << "/" << table.size() << ", plain_sum=" << workspace->sum_plain;
      gsl_sum_levin_utrunc_free(workspace);
   }
   return sum;
}

// --------------------------------------------------------------------------------------------------------------------------------

double root_fsolver_wrapper::findRoot(double tol_abs, double tol_rel, char const* funcName, uint maxIter) const
{
   auto s = solver_.get();
   for (uint i = 0;; ++i)
   {
      gsl_root_fsolver_iterate(s);
      double low = gsl_root_fsolver_x_lower(s);
      double high = gsl_root_fsolver_x_upper(s);
      int status = gsl_root_test_interval(low, high, tol_abs, tol_rel);
      if (status != GSL_CONTINUE) break;
      if (i >= maxIter) throw std::runtime_error(std::string("findRoot failed to converge in ") + std::string(funcName));
   }
   return gsl_root_fsolver_root(s);
}

// --------------------------------------------------------------------------------------------------------------------------------


#if defined(ENABLE_GF_TESTFUNCTIONS)

void GF_EXPORT TestSphericalBessel(uint n, uint size, double zstart, double dz, double* j_buffer, double* y_buffer)
{
   const SphericalBesselGenerator& sbt(SphericalBesselGenerator::instance());
   for (uint i = 0; i < size; ++i)
   {
      double z = zstart + i*dz;
      j_buffer[i] = sbt.j(n, z);
      y_buffer[i] = sbt.y(n, z);
   }
}

void GF_EXPORT TestCylindricalBessel(uint n, uint size, double zstart, double dz, double* J_buffer, double* Y_buffer)
{
   const CylindricalBesselGenerator& cbt(CylindricalBesselGenerator::instance());
   for (uint i = 0; i < size; ++i)
   {
      double z = zstart + i*dz;
      J_buffer[i] = cbt.J(n, z);
      Y_buffer[i] = cbt.Y(n, z);
   }
}

#include "GreensFunction.hpp"
#include "PairGreensFunction.hpp"
#include "GreensFunction3D.hpp"
#include "GreensFunction3DAbs.hpp"
#include "GreensFunction3DAbsSym.hpp"
#include "GreensFunction3DRadAbs.hpp"
#include "GreensFunction3DRadInf.hpp"
#include "math.h"

void GF_EXPORT TestGreensFunction(uint gfn, uint gfc, uint size, double *buffer, double D, double r0, double a, double kf, double sigma, double t, double r)
{
   gsl_set_error_handler([](const char*, const char*, int, int) { throw  std::runtime_error("GSL"); });

   std::unique_ptr<GreensFunction> gf1;
   std::unique_ptr<PairGreensFunction> gf2;
   switch (gfn)
   {
   case 0: gf1 = std::make_unique<GreensFunction3DAbsSym>(GreensFunction3DAbsSym(D, a)); break;
   case 1: gf2 = std::make_unique<GreensFunction3D>(GreensFunction3D(D, r0)); break;
   case 2: gf2 = std::make_unique<GreensFunction3DAbs>(GreensFunction3DAbs(D, r0, a)); break;
   case 3: gf2 = std::make_unique<GreensFunction3DRadAbs>(GreensFunction3DRadAbs(D, kf, r0, sigma, a)); break;
   case 4: gf2 = std::make_unique<GreensFunction3DRadInf>(GreensFunction3DRadInf(D, kf, r0, sigma)); break;
   }

   if (gfc < 3 && gfn > 0) gf1 = std::move(gf2);
   if (gfc == 3 && gf2.get() == nullptr) return;

   for (uint i = 0; i < size; ++i)
   {
      double rnd = (double)i / size;

      try
      {
         switch (gfc)
         {
         case 0: buffer[i] = gf1->drawTime(rnd); break;
         case 1: buffer[i] = gf1->drawR(rnd, t); break;
         case 2: buffer[i] = (double)(int)gf1->drawEventType(rnd, t); break;
         case 3: buffer[i] = gf2->drawTheta(rnd, r, t); break;
         }
      }
      catch (...)
      {
         buffer[i] = NAN;
      }
   }
}

double GF_EXPORT DrawGreensFunction(uint gfn, uint gfc, double rnd, double D, double r0, double a, double kf, double sigma, double t, double r)
{
   std::unique_ptr<GreensFunction> gf1;
   std::unique_ptr<PairGreensFunction> gf2;
   switch (gfn)
   {
   case 0: gf1 = std::make_unique<GreensFunction3DAbsSym>(GreensFunction3DAbsSym(D, a)); break;
   case 1: gf2 = std::make_unique<GreensFunction3D>(GreensFunction3D(D, r0)); break;
   case 2: gf2 = std::make_unique<GreensFunction3DAbs>(GreensFunction3DAbs(D, r0, a)); break;
   case 3: gf2 = std::make_unique<GreensFunction3DRadAbs>(GreensFunction3DRadAbs(D, kf, r0, sigma, a)); break;
   case 4: gf2 = std::make_unique<GreensFunction3DRadInf>(GreensFunction3DRadInf(D, kf, r0, sigma)); break;
   }

   if (gfc < 3 && gfn > 0) gf1 = std::move(gf2);
   if (gfc == 3 && gf2.get() == nullptr) return NAN;

   try
   {
      switch (gfc)
      {
      case 0: return gf1->drawTime(rnd);
      case 1: return gf1->drawR(rnd, t);
      case 2: return (double)(int)gf1->drawEventType(rnd, t);
      case 3: return gf2->drawTheta(rnd, r, t);
      }
   }
   catch (...)
   {
   }
   return NAN;
}


#endif

// --------------------------------------------------------------------------------------------------------------------------------
