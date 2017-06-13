#ifndef GREENSFUNCTION2DRADABS_HPP
#define GREENSFUNCTION2DRADABS_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <array>
#include "PairGreensFunction.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// Greens function class for 2d Green's Function for 2d annulus with radial and
// axial dependence. Inner boundary is radiative (rad) (reaction event), outer 
// boundary is absorbing (abs) (escape event). Different "draw" functions 
// provide a way to draw certain values from the Green's Function, e.g. an
// escape angle theta ("drawTheta" function).
// 
// Written by Laurens Bossen. Adapted by Martijn Wehrens.
// FOM Institute AMOLF.

class GF_EXPORT GreensFunction2DRadAbs : public PairGreensFunction
{
   static const uint MAX_ALPHA_SEQ;            // The maximum number of n terms
   static const double SCAN_START;               // Left boundary of 1st search interval 1st root
   static const double FRACTION_SCAN_INTERVAL;   // Length of the scanning interval relative to estimated interval TODO CHANGED THIS FROM .5 to .2
   static const double CONVERGENCE_ASSUMED;
   static const double INTERVAL_MARGIN;

   // After CONVERGENCE_ASSUMED subsequent roots that lay within +/- 
   // INTERVAL_MARGIN from the distance to which the distance is known to 
   // converge, it is assumed all following roots have a distances in-between
   // that don't deviate for more than INTERVAL_MARGIN from the distance to 
   // which the roots are known to converge (Pi/(a-sigma)).

public:

   GreensFunction2DRadAbs(double D, double kf, double r0, double sigma, double a);

   virtual std::string dump() const override;

   virtual const char* type_name() const override { return "GreensFunction2DRadAbs"; }

   virtual double drawR(double rnd, double t) const override;

   virtual double drawTheta(double rnd, double r, double t) const override;

   double geth() const { return h_; }

   double geta() const { return a_; }

   double drawTime(double rnd) const override;

   EventKind drawEventType(double rnd, double t) const override;

   double f_alpha0(double alpha) const;

   double f_alpha(double alpha, int n) const;

   double p_survival(double t) const;

   double p_survival_table(double t, DoubleVector& table) const;

   double leaves(double t) const;

   double leavea(double t) const;

   double p_m(int n, double r, double t) const;

   double dp_m_at_a(int m, double t) const;

   double p_m_alpha(const uint n, uint m, double r, double t) const;

   double dp_m_alpha_at_a(const uint n, uint m, double t) const;

   void GiveRootInterval(double& low, double& high, int n) const;

   void GiveRootIntervalSimple(double& low, double& high, int n, uint i) const;

   double getAlphaRoot0(double low, double high) const;

   double getAlphaRootN(double low, double high, int n) const;

   double getAlphaRoot(double high, double low, int n) const;

   void decideOnMethod2(uint n, uint i) const;

   double getAlpha(uint n, uint i) const;

   double p_survival_i(double alpha) const;

   double calc_A_i_0(double alpha) const;

   double leaves_i(double alpha) const;

   std::tuple<double, double, double> Y0J0J1_constants(double alpha, double t) const;

   double givePDFTheta(double theta, double r, double t) const;

   double givePDFR(double r, double t) const;

   void dumpRoots(int n) const;

private:

   void clearAlphaTable() const;

   DoubleVector& getAlphaTable(uint n) const
   {
      THROW_UNLESS(std::invalid_argument, n < alphaTables_.size());
      return alphaTables_[n];
   }

   double p_int_r_table(double r, const DoubleVector& Y0_aAnTable, const DoubleVector& J0_aAnTable, const DoubleVector& Y0J1J0Y1Table) const;

   double ip_theta_table(double theta, const DoubleVector& p_nTable) const;

   double p_survival_i_exp_table(uint i, double t, const DoubleVector& table) const;

   double leavea_i_exp(uint i, double alpha) const;

   double leaves_i_exp(uint i, double alpha) const;

   double ip_theta_n(uint m, double theta, const DoubleVector& p_nTable) const;

   double p_int_r_i_exp_table(uint i, double r, const DoubleVector& Y0_aAnTable, const DoubleVector& J0_aAnTable, const DoubleVector& Y0J1J0Y1Table) const;

   void createPsurvTable(DoubleVector& table) const;

   void createY0J0Tables(DoubleVector& Y0_Table, DoubleVector& J0_Table, DoubleVector& Y0J1J0Y1_Table, double t) const;

   void makep_mTable(DoubleVector& p_mTable, double r, double t) const;

   void makedp_m_at_aTable(DoubleVector& p_mTable, double t) const;

   uint guess_maxi(double t) const;

   struct p_survival_table_params
   {
      const GreensFunction2DRadAbs& gf;
      DoubleVector& table;
      double rnd;
   };

   static double p_survival_table_F(double t, const p_survival_table_params* const params);

   struct p_int_r_params
   {
      const GreensFunction2DRadAbs& gf;
      const DoubleVector& Y0_aAnTable;
      const DoubleVector& J0_aAnTable;
      const DoubleVector& Y0J1J0Y1Table;
      double rnd;
   };

   static double p_int_r_F(double r, const p_int_r_params* const params);


   double ip_theta_F(double theta, const DoubleVector& p_nTable, double value) const;

   const double h_;
   const double a_;

   // Tables that hold calculated roots (y=0) of "alpha" function for each
   // order n.
   mutable std::array<DoubleVector, GfCfg.MAX_ORDER()> alphaTables_;

   // Constants used in the roots of f_alpha() finding algorithm.
   // ====
   //
   // This constant will simply be M_PI/(a-Sigma), the value to which the 
   // distance between roots of f_alpha() should converge.
   double estimated_alpha_root_distance_;
   //
   // Table which tells us at which x we're left with scanning the alpha 
   // function for a sign change, for a given order n. (A sign change would 
   // indicate a root (y=0) lies between the boundaries of the "scanned" 
   // interval.) 
   //      If x_scan[n] < 0, this indicates scanning is no longer required
   // because the distance between the roots is converging and within 
   // boundaries that allow the direct use of the estimate interval width 
   // pi/(sigma-a).
   //      Initial values are set by constructor.
   mutable std::array<double, GfCfg.MAX_ORDER()> alpha_x_scan_table_;
   //
   // Table that keeps track of the number of previous subsequent roots that 
   // we're within margin of the distance to which they're expected to 
   // converge.
   mutable std::array<int, GfCfg.MAX_ORDER()> alpha_correctly_estimated_;

};

// --------------------------------------------------------------------------------------------------------------------------------

#endif // GREENSFUNCTION2DRADABS_HPP