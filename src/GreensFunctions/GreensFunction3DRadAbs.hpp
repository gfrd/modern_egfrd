#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <array>
#include <functional>
#include "helperFunctionsGf.hpp"
#include "PairGreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GF_EXPORT GreensFunction3DRadAbs : public PairGreensFunction
{
   static const uint MAX_ALPHA_SEQ;
   static const double H;

public:

   GreensFunction3DRadAbs(double D, double kf, double r0, double sigma, double a);

   std::string dump() const override;

   const char* type_name() const override { return "GreensFunction3DRadAbs"; }

   double geta() const { return a_; }

   double geth() const { return h_; }

   double drawTime(double rnd) const override;

   double drawR(double rnd, double t) const override;

   EventKind drawEventType(double rnd, double t) const override;

   double drawTheta(double rnd, double r, double t) const override;


private:

   friend class UnitTestGreensFunction3DRadAbs; // for UnitTests

   double alpha0_i(int i) const;
   
   double f_alpha0(double alpha) const;

   double ip_theta(double theta, double r, double t) const;

   double f_alpha0_aux(double alpha) const;
   double f_alpha(double alpha, int n) const;
   double f_alpha_aux(double alpha, int n) const;
   double p_survival_table(double t, DoubleVector& table) const;
   double leaves(double t) const;
   double leavea(double t) const;
   double p_int_r(double r, double t) const;
   double p_n_alpha(uint i, uint n, double r, double t) const;
   double dp_n_alpha_at_a(uint i, uint n, double t) const;
   uint alphaOffset(uint n) const;
   double alpha_i(int i, int n, const root_fsolver_wrapper& solver) const;
   double p_survival_i(double alpha) const;
   double leavea_i(double alpha) const;
   double leaves_i(double alpha) const;
   double p_int_r_i(double r, double alpha, double num_r0) const;
   void clearAlphaTable() const;
   DoubleVector& getAlphaTable(uint n) const { THROW_UNLESS(std::invalid_argument, n < alphaTables_.size()); return alphaTables_[n]; }
   double getAlpha(uint n, uint i) const;
   double getAlpha0(uint i) const;
   double p_survival_i_exp_table(uint i, double t, const DoubleVector& table) const;
   double leavea_i_exp(uint i, double alpha) const;
   double leaves_i_exp(uint i, double alpha) const;
   double p_int_r_i_exp(uint i, double t, double r) const;
   void updateAlphaTable0(double t) const;
   void updateAlphaTable(uint n, double t) const;
   void createPsurvTable(DoubleVector& table) const;
   uint guess_maxi(double t) const;
   double num_r0(double alpha) const;
   double funcSumMaxAlpha(int n, double max_alpha, const std::function<double(uint, uint)>& func) const;
   void make_pn_table(DoubleVector& p_nTable, double t, double factor, const std::function<double(uint, uint)>& func, bool max_Dt = false) const;

   const double a_;
   const double h_;
   const double hsigmap1_;

   mutable std::array<int, GfCfg.MAX_ORDER()+1> alphaOffsetTables_;
   mutable std::array<DoubleVector, GfCfg.MAX_ORDER()+1> alphaTables_;
};

// --------------------------------------------------------------------------------------------------------------------------------
