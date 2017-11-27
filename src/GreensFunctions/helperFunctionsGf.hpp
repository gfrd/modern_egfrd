#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <functional>
#include <memory>
#include <gsl/gsl_roots.h>
#include "DefsGf.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

double funcSum_all(std::function<double(uint i)> f, uint max_i);

double funcSum(std::function<double(uint i)> f, uint max_i, double tolerance = GfCfg.TOLERANCE);

// --------------------------------------------------------------------------------------------------------------------------------

std::string GetExecutablePath();

bool file_exists(const std::string& name);

// --------------------------------------------------------------------------------------------------------------------------------

// Template class to map Lambda functions to GSL function pointers
// Uses params as this pointer, so no params avail to functions
template<typename TLambda>
class gsl_lambda : public gsl_function
{
public:
   gsl_lambda(const TLambda& func) : func_(func)
   {
      function = &gsl_lambda::invoke;
      params = this;
   }

private:
   const TLambda& func_;
   static double invoke(double x, void *params)
   {
      return static_cast<gsl_lambda*>(params)->func_(x);
   }
};

// --------------------------------------------------------------------------------------------------------------------------------

class root_fsolver_wrapper
{
public:
   root_fsolver_wrapper(const gsl_root_fsolver_type* type = gsl_root_fsolver_brent) : solver_(gsl_root_fsolver_alloc(type), gsl_root_fsolver_free) { }

   void set(gsl_function *f, double x_lower, double x_upper) const
   {
      gsl_root_fsolver_set(solver_.get(), f, x_lower, x_upper);
   }

   double findRoot(double tol_abs, double tol_rel, char const* funcName, uint maxIter = 100) const;

private:
   std::unique_ptr<gsl_root_fsolver, void(*)(gsl_root_fsolver*)> solver_;
};

// --------------------------------------------------------------------------------------------------------------------------------

//#define ENABLE_GF_TESTFUNCTIONS
#if defined(ENABLE_GF_TESTFUNCTIONS)
extern "C"
{
   void GF_EXPORT TestSphericalBessel(uint n, uint size, double z, double dz, double* j_buffer, double* y_buffer);
   void GF_EXPORT TestCylindricalBessel(uint n, uint size, double z, double dz, double* J_buffer, double* Y_buffer);
   void GF_EXPORT TestGreensFunction(uint gfn, uint gfc, uint size, double *buffer, double D, double r0, double a, double kf, double sigma, double t, double r);
   double GF_EXPORT DrawGreensFunction(uint gfn, uint gfc, double rnd, double D, double r0, double a, double kf, double sigma, double t, double r);
}
#endif

// --------------------------------------------------------------------------------------------------------------------------------
