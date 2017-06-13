#ifndef GREENSFUNCTION3DRADINF_HPP
#define GREENSFUNCTION3DRADINF_HPP 

// --------------------------------------------------------------------------------------------------------------------------------

#include "PairGreensFunction.hpp"
#include <gsl/gsl_integration.h>

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction3DRadInf : public PairGreensFunction
{
   struct p_corr_R_params;
   struct p_theta_params;

public:

   GreensFunction3DRadInf(double D, double kf, double r0, double sigma);

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction3DRadInf"; }

   double getkD() const { return kD_; }

   double getalpha() const { return alpha_; }

   double drawTime(double rnd) const override;

   double drawR(double rnd, double t) const override;

   double drawTheta(double rnd, double r, double t) const override;

private:

   friend class UnitTestGreensFunction3DRadInf;  // for UnitTests

   double p_int_r(double r, double t) const;

   double p_corr_R(double alpha, uint n, double r, double t) const;

   double ip_corr_table(double theta, double r, const DoubleVector& RnTable) const;

   double ip_theta_table(double r, double theta, double time, const DoubleVector& RnTable) const;

   DoubleVector makeRnTable(double r, double t) const;

   double Rn(uint order, double r, double t, gsl_integration_workspace* workspace, double tol) const;

   const double kD_;
   const double alpha_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif // GREENSFUNCTION3DRADINF_HPP