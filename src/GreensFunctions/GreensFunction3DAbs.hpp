#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "PairGreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction3DAbs : public PairGreensFunction
{
   static const double MIN_T;
   static const uint MAX_ALPHA_SEQ;

public:

   GreensFunction3DAbs(double D, double r0, double a);

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction3DAbs"; }

   double geta() const { return a_; }

   double drawTime(double rnd) const override;

   double drawR(double rnd, double t) const override;

   double drawTheta(double rnd, double r, double t) const override;

private:

   friend class UnitTestGreenFunction3DAbs; // for UnitTests

   double p_survival(double t) const;

   double dp_survival(double t) const;

   double p_int_r(double r, double t) const;

   double p_theta(double theta, double r, double t) const;

   double ip_theta(double theta, double r, double t) const;

   double dp_theta(double theta, double r, double t) const;

   double idp_theta(double theta, double r, double t) const;

   double p_n(int n, double r, double t) const;

   double dp_n(int n, double t) const;

   double p_n_alpha(uint i, uint n, double r, double t) const;

   double dp_n_alpha(uint i, uint n, double t) const;

   double p_theta_table(double theta, const DoubleVector& p_nTable) const;

   double ip_theta_table(double theta, const DoubleVector& p_nTable) const;

   void makep_nTable(DoubleVector& p_nTable, double r, double t) const;

   void makedp_nTable(DoubleVector& p_nTable, double t) const;

   const double a_;
};

// --------------------------------------------------------------------------------------------------------------------------------
