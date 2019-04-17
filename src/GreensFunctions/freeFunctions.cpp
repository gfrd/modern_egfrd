#include <cmath>
#include <algorithm>
#include <stdexcept>
#include <gsl/gsl_errno.h>
#include "freeFunctions.hpp"
#include "helperFunctionsGf.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

//   Calculates std::exp(x^2) * std::erfc(x)
//   See asymptotic expansion here:
//   http://en.wikipedia.org/wiki/Error_function

double expxsq_erfc(double x)
{
   double xsq = x * x;
   if (x > 26.0)
   {
      double M_1_SQRTPI = M_2_SQRTPI * 0.5;
      double x2sq_r = 1.0 / (2.0 * xsq);  // 2 / (2 x)^2
      return M_1_SQRTPI / x * (1.0 - x2sq_r + x2sq_r * x2sq_r);
   }

   return std::exp(xsq) * std::erfc(x);
}

// --------------------------------------------------------------------------------------------------------------------------------

// W(a, b) = std::exp(2 a b + b^2) std::erfc(a + b)
double W(double a, double b)
{
   // std::exp(2 a b + b^2) std::erfc(a + b) == std::exp(- a^2) std::exp((a + b)^2) std::erfc(a + b)
   return std::exp(-a * a) * expxsq_erfc(a + b);
}

// --------------------------------------------------------------------------------------------------------------------------------

/* List of 1D functions used as approximations for the greensfunctions
   of the full domain. For small t, a particle only scans part of the domain
   and thus not all the boundary conditions have to be included. Because the
   greensfunctions converge bad for small t, we use these approximations.

   note: drift not included yet!

   *** Naming:
   All 1D function names start with an X, followed by:

   P for position density function.
   S for survival probability function.
   I for cumulative distribution function of position.

   ** Last two number give boundaries
   (By J.V. Beck - Heat Conduction in Green's Functions)

   0 - No boundary.
   1 - Absorbing boundary.
   2 - Reflective boundary.
   3 - Radiation boundary with rate ka.

   The approximation for the greensfunction with a sink has 3 boundaries: 030
   left absorbing boundary, sink, right absorbing boundary.
   */

   // --------------------------------------------------------------------------------------------------------------------------------

double XP00(double r, double t, double r0, double D, double v)
{
   double fourDt = 4 * D * t;
   double rminr0minv2 = gsl_pow_2(r - r0 - v * t);

   return 1.0 / std::sqrt(fourDt * M_PI) * std::exp(-rminr0minv2 / fourDt);
}

// --------------------------------------------------------------------------------------------------------------------------------

double XI00(double r, double t, double r0, double D, double v)
{
   double sqrt4Dt = std::sqrt(4 * D * t);
   double rminr0minv = r - r0 - v * t;
   return 0.5 * (1.0 + std::erf(rminr0minv / sqrt4Dt));
}

// --------------------------------------------------------------------------------------------------------------------------------

double XS00(double, double, double, double)
{
   return 1.0;
}

// --------------------------------------------------------------------------------------------------------------------------------

double XP10(double r, double t, double r0, double D, double v)
{
   double fourDt = 4 * D * t;
   double rminr02 = gsl_pow_2(r - r0);
   double rplusr02 = gsl_pow_2(r + r0);
   double v_2 = v / 2.0;
   double drift_prefac = v != 0.0 ? 1.0 : std::exp(v_2 / D * (r - r0 - v_2 * t));
   return drift_prefac / std::sqrt(fourDt * M_PI) * (std::exp(-rminr02 / fourDt) - std::exp(-rplusr02 / fourDt));
}

// --------------------------------------------------------------------------------------------------------------------------------

double XI10(double r, double t, double r0, double D, double v)
{
   double sqrt4Dt = std::sqrt(4 * D * t);
   if (v == 0.0)
   {
      double rminr0 = r - r0;
      double rplusr0 = r + r0;
      return 0.5 * (2.0 * std::erf(r0 / sqrt4Dt) + std::erf(rminr0 / sqrt4Dt) - std::erf(rplusr0 / sqrt4Dt));
   }

   double r0minvt = r0 - v*t;
   double r0plusvt = r0 + v*t;
   double term1 = std::erfc((r0minvt + r) / sqrt4Dt) - std::erfc(r0minvt / sqrt4Dt);
   double term2 = std::erf(r0plusvt / sqrt4Dt) - std::erf((r0plusvt - r) / sqrt4Dt);
   return 0.5 * (std::exp(-v * r0 / D) * term1 + term2);
}

// --------------------------------------------------------------------------------------------------------------------------------

double XS10(double t, double r0, double D, double v)
{
   double sqrt4Dt = std::sqrt(4 * D * t);
   if (v == 0.0) return std::erf(r0 / sqrt4Dt);

   double vt = v * t;
   double corr = std::exp((6 * v*r0 - v*vt) / (4 * D));           //std::erfc( -(r0 + vt) / sqrt4Dt ) - std::exp(- v * r0 / D ) * std::erfc( (r0 - vt) / sqrt4Dt) );
   return 0.5 * corr * W(r0 / sqrt4Dt, -vt / sqrt4Dt);
}

// --------------------------------------------------------------------------------------------------------------------------------

double XP20(double r, double t, double r0, double D, double v)
{
   double fourDt = 4 * D * t;
   double rminr02 = gsl_pow_2(r - r0);
   double rplusr02 = gsl_pow_2(r + r0);
   double XP20_nov = 1.0 / std::sqrt(fourDt * M_PI) * (std::exp(-rminr02 / fourDt) + std::exp(-rplusr02 / fourDt));
   if (v == 0.0) return XP20_nov;

   double sqrt4Dt = fourDt;
   double v_2 = v / 2.0;
   double drift_prefac = std::exp(v_2 / D * (r - r0 - v_2 * t));
   double XP20_v = std::exp(v_2 / D * (r + r0) + v_2 * v_2 / D * t) * std::erfc((r + r0 + v * t) / sqrt4Dt);
   return drift_prefac * (v_2 / D * XP20_v + XP20_nov);
}

// --------------------------------------------------------------------------------------------------------------------------------

double XI20(double r, double t, double r0, double D, double v)
{
   double sqrt4Dt = std::sqrt(4 * D * t);
   if (v == 0.0) return 0.5 * std::erf((r - r0) / sqrt4Dt) + std::erf((r + r0) / sqrt4Dt);

   throw std::runtime_error("TODO: c.d.f with drift not included yet. Does exist!");
}

// --------------------------------------------------------------------------------------------------------------------------------

double XS20(double t, double r0, double D, double v)
{
   UNUSED(t);
   UNUSED(r0);
   UNUSED(D);
   UNUSED(v);
   return 1.0;
}

// --------------------------------------------------------------------------------------------------------------------------------

double XP30term_nov(double r, double t, double r0, double ka, double D)
{
   double r0a = std::fabs(r0);
   double fourDt = 4 * D * t;
   double k_D = ka / D;
   double rplusr02 = gsl_pow_2(r + r0a);
   double arg = (r + r0a) / std::sqrt(fourDt) + ka * std::sqrt(t / D);
   return -k_D * std::exp(-rplusr02 / fourDt) * expxsq_erfc(arg);
}

// --------------------------------------------------------------------------------------------------------------------------------

double XP30term_v(double r, double t, double r0, double ka, double D, double v)
{
   double r0a = std::fabs(r0);
   double sqrt4Dt = std::sqrt(4 * D * t);
   double v_2 = v / 2.0;
   double kplusv2 = ka + v_2;
   double erfc_arg = (r + r0a + 2 * kplusv2 * t) / sqrt4Dt;
   double exp_arg = 1.0 / D * (kplusv2 * kplusv2 * t + kplusv2 * (r + r0a));
   return -std::exp(exp_arg) * std::erfc(erfc_arg);
}

// --------------------------------------------------------------------------------------------------------------------------------

double XP30(double r, double t, double r0, double ka, double D, double v)
{
   double r0a = std::fabs(r0);
   double XP20temp = XP20(r, t, r0a, D, 0.0);
   if (v == 0.0) return XP20temp + XP30term_nov(r, t, r0a, ka, D);

   double v_2 = v / 2.0;
   double drift_prefac = std::exp(v_2 / D * (r - r0a - v_2 * t));
   return drift_prefac * (XP20temp + 1 / D * (ka + v / 2) * XP30term_v(r, t, r0a, ka, D, v));
}

// --------------------------------------------------------------------------------------------------------------------------------

double XI30term_nov(double r, double t, double r0, double ka, double D)
{
   double sqrt4Dt = std::sqrt(4 * D * t);
   double r0a = std::fabs(r0);

   double term1 = std::erf(r0a / sqrt4Dt) - std::erf((r + r0a) / sqrt4Dt);
   double term2 = W(r0a / sqrt4Dt, 2 * ka*t / sqrt4Dt);           //std::exp( k_D * (ka * t + r0a) ) * std::erfc( (2 * ka * t + r0a)/sqrt4Dt ) );
   double term3 = W((r + r0a) / sqrt4Dt, 2 * ka*t / sqrt4Dt);     //std::exp( k_D * (ka * t + r0a + r) ) * std::erfc( (2 * ka * t + r0a + r)/sqrt4Dt ) );
   return term1 + term2 - term3;
}

// --------------------------------------------------------------------------------------------------------------------------------

double XI30(double r, double t, double r0, double ka, double D, double v)
{
   double r0a = std::fabs(r0);
   if (v == 0.0) return XI20(r, t, r0a, D, 0.0) + XI30term_nov(r, t, r0a, ka, D);

   throw std::runtime_error("TODO: c.d.f with drift not included yet. Does exist!");
}

// --------------------------------------------------------------------------------------------------------------------------------

double XS30(double t, double r0, double ka, double D, double v)
{
   double r0a = std::fabs(r0);
   double sqrt4Dt = std::sqrt(4 * D * t);
   if (v == 0.0) return std::erf(r0a / sqrt4Dt) + W(r0a / sqrt4Dt, 2 * ka * t / sqrt4Dt);    // std::exp( k_D * ka * t + k_D * r0a ) * std::erfc( (2 * ka * t + r0a) / sqrt4Dt );

   double k_D = ka / D;
   double v_2 = v / 2.0;
   double kplusv2 = ka + v_2;
   double r0plusvt = r0a + v * t;
   double r0minvt = r0a - v * t;
   double erfc_arg = (r0a + 2 * kplusv2 * t) / sqrt4Dt;
   double exp_arg = k_D * (r0a + (ka + v) * t);
   double term2 = std::erfc(-r0plusvt / sqrt4Dt) - ka / (ka + v) * std::exp(-v / D * r0a) * std::erfc(r0minvt / sqrt4Dt);
   return (ka + v_2) / (ka + v) * std::exp(exp_arg) * std::erfc(erfc_arg) + 0.5 * term2;
}

// --------------------------------------------------------------------------------------------------------------------------------

double XP030(double r, double t, double r0, double ka, double D)
{
   double r0a = std::fabs(r0);
   return XP00(r, t, r0a, D, 0.0) + 0.5 * XP30term_nov(std::fabs(r), t, r0a, 0.5 * ka, D);
}

// --------------------------------------------------------------------------------------------------------------------------------

double XI030(double r, double t, double r0, double ka, double D)
{
   double ka_2 = ka * 0.5;
   double r0a = std::fabs(r0);
   double sign = r < 0 ? -1 : 1;
   return XI00(r, t, r0a, D, 0.0) + 0.5 * ((XS30(t, r0a, ka_2, D, 0.0) - 1) + sign * XI30term_nov(std::fabs(r), t, r0a, ka_2, D));
}

// --------------------------------------------------------------------------------------------------------------------------------

double XS030(double t, double r0, double ka, double D)
{
   return XS30(t, std::fabs(r0), 0.5 * ka, D, 0.0);
}

// --------------------------------------------------------------------------------------------------------------------------------

double __p_irr(double r, double t, double r0, double D, double sigma, double alpha)
{
   double sqrtD = std::sqrt(D);
   double Dt4 = 4.0 * D * t;
   double r_plus_r0_minus_2sigma = r + r0 - 2.0 * sigma;
   double num1 = std::exp(-gsl_pow_2(r - r0) / Dt4);
   double num2 = std::exp(-gsl_pow_2(r_plus_r0_minus_2sigma) / Dt4);
   double num3 = W(r_plus_r0_minus_2sigma / std::sqrt(Dt4), alpha * std::sqrt(t));
   double num = (num1 + num2) / std::sqrt(4.0 * M_PI * t) - alpha * num3;
   double den = 4.0 * M_PI * r * r0 * sqrtD;
   return 4.0 * M_PI * r * r * num / den;
}

// --------------------------------------------------------------------------------------------------------------------------------

double p_irr(double r, double t, double r0, double kf, double D, double sigma)
{
   double kD = 4.0 * M_PI * sigma * D;
   double alpha = (1.0 + (kf / kD)) * (std::sqrt(D) / sigma);
   return __p_irr(r, t, r0, D, sigma, alpha);
}

// --------------------------------------------------------------------------------------------------------------------------------

double p_survival_irr(double t, double r0, double kf, double D, double sigma)
{
   double kD = 4.0 * M_PI * sigma * D;
   double alpha = (1.0 + kf / kD) * (std::sqrt(D) / sigma);
   double p = __p_reaction_irr(t, r0, kf, D, sigma, alpha, kD);
   return 1.0 - p;
}

// --------------------------------------------------------------------------------------------------------------------------------

double __p_reaction_irr(double t, double r0, double kf, double D, double sigma, double alpha, double kD)
{
   double sqrtt = std::sqrt(t);
   double sqrtD = std::sqrt(D);
   double r0_m_sigma_over_sqrt4D_t = (r0 - sigma) / ((sqrtD + sqrtD) * sqrtt);
   double Wf = W(r0_m_sigma_over_sqrt4D_t, alpha * sqrtt);
   double factor = sigma * kf / (r0 * (kf + kD));
   return factor * (std::erfc(r0_m_sigma_over_sqrt4D_t) - Wf);
}

// --------------------------------------------------------------------------------------------------------------------------------

double __p_reaction_irr_t_inf(double r0, double kf, double sigma, double kD)
{
   double kf_kD_r0 = (kf + kD) * r0;
   return 1 - (kf_kD_r0 - kf * sigma) / kf_kD_r0;
}

// --------------------------------------------------------------------------------------------------------------------------------

double p_survival_nocollision(double t, double r0, double D, double a)
{
   double Dt = D * t;
   double asq = a * a;
   double a_r = 1.0 / a;
   double asq_r = a_r * a_r;
   double PIr0 = M_PI * r0;

   double angle_factor = PIr0 * a_r;
   double exp_factor = -Dt * M_PI * M_PI * asq_r;

   uint i_max = std::max(static_cast<uint>(std::ceil(std::sqrt(M_PI * M_PI + asq * std::log(1.0 / GfCfg.TOLERANCE) / Dt) * M_1_PI)), 2u);

   double p = 0.0;
   double sign = 1.0;
   for (uint i = 1; i <= i_max; ++i)
   {
      p += sign * std::exp(exp_factor * i * i) * std::sin(angle_factor * i) / i;
      sign = -sign;
   }
   return p * ((a + a) / PIr0);
}

// --------------------------------------------------------------------------------------------------------------------------------

double dp_survival_nocollision(double t, double r0, double D, double a)
{
   double Dt = D * t;
   double asq = a * a;
   double a_r = 1.0 / a;
   double asq_r = a_r * a_r;
   double PIr0 = M_PI * r0;

   double angle_factor = PIr0 * a_r;
   double exp_factor = -Dt * M_PI * M_PI * asq_r;

   uint i_max = std::max(static_cast<uint>(std::ceil(std::sqrt(M_PI * M_PI + asq * std::log(1.0 / GfCfg.TOLERANCE) / Dt) * M_1_PI)), 2u);

   double p = 0.0;
   double sign = -1.0;
   for (uint i = 1; i <= i_max; ++i)
   {
      p += sign * std::exp(exp_factor * i * i) * std::sin(angle_factor * i) * i;
      sign = -sign;
   }

   double factor = D * (M_PI + M_PI) / (a * r0);
   return p * factor;
}

// --------------------------------------------------------------------------------------------------------------------------------

double p_theta_free(double theta, double r, double r0, double t, double D)
{
   double Dt4 = 4.0 * D * t;
   double Dt4Pi = Dt4 * M_PI;
   double term1 = std::exp(-(r * r - 2.0 * std::cos(theta) * r * r0 + r0 * r0) / Dt4);
   double term2 = 1.0 / std::sqrt(Dt4Pi * Dt4Pi * Dt4Pi);
   return term1 * term2 * std::sin(theta);
}

// --------------------------------------------------------------------------------------------------------------------------------

double ip_theta_free(double theta, double r, double r0, double t, double D)
{
   double Dt = D * t;
   double Dt2 = Dt + Dt;
   double rr0 = r * r0;
   double rr0_over_2Dt = rr0 / Dt2;
   double rsqr0sq_over_4Dt = (r * r + r0 * r0) / (Dt2 + Dt2);

   double term1 = std::expm1(rr0_over_2Dt - rsqr0sq_over_4Dt);
   double term2 = std::expm1(rr0_over_2Dt * std::cos(theta) - rsqr0sq_over_4Dt);

   double den = 4.0 * std::sqrt(M_PI * M_PI * M_PI * Dt) * rr0;
   return (term1 - term2) / den;
}

// --------------------------------------------------------------------------------------------------------------------------------

double g_bd_3D(double r, double sigma, double t, double D)
{
   double Dt4 = 4.0 * D * t;
   double mDt4_r = -1.0 / Dt4;
   double sqrtDt4 = std::sqrt(Dt4);
   double sqrtDt4_r = 1.0 / sqrtDt4;

   double rps = r + sigma;
   double rms = r - sigma;

   double term1 = (std::exp(rps * rps * mDt4_r) - std::exp(rms * rms * mDt4_r)) * sqrtDt4 / (M_SQRTPI * r);
   double term2 = std::erf(rps * sqrtDt4_r) - std::erf(rms * sqrtDt4_r);
   return 0.5 * (term1 + term2) * r * r;
}

// --------------------------------------------------------------------------------------------------------------------------------

double g_bd_1D(double r, double sigma, double t, double D, double v)
{
   double Dt4 = 4.0 * D * t;
   double sqrtDt4 = std::sqrt(Dt4);
   double sqrtDt4_r = 1.0 / sqrtDt4;
   const double vt = v*t;

   double s_plus_r_plus_vt = sigma + r + vt;
   double s_min_r_min_vt = sigma - r - vt;

   return 0.5 * std::erf(s_plus_r_plus_vt * sqrtDt4_r) + std::erf(s_min_r_min_vt * sqrtDt4_r);
}

// --------------------------------------------------------------------------------------------------------------------------------

double I_bd_3D(double sigma, double t, double D)
{
   double Dt = D * t;
   double Dt2 = Dt + Dt;
   double sqrtDt = std::sqrt(Dt);
   double sigmasq = sigma * sigma;

   double term1 = 1.0 / (3.0 * M_SQRTPI);
   double term2 = sigmasq - Dt2;
   double term3 = Dt2 - 3.0 * sigmasq;
   double term4 = M_SQRTPI * sigmasq * sigma * std::erfc(sigma / sqrtDt);

   return term1 * (-sqrtDt * (term2 * std::exp(-sigmasq / Dt) + term3) + term4);
}

// --------------------------------------------------------------------------------------------------------------------------------

double I_bd_1D(double sigma, double t, double D, double v)
{
   if (D == 0) return 0;

   double Dt4 = 4 * D * t;
   double sqrt4Dt = std::sqrt(Dt4);
   double vt = v*t;

   double arg1 = -(2 * sigma + vt)*(2 * sigma + vt) / Dt4;
   double term1 = std::exp(-vt*vt / Dt4) - std::exp(arg1);
   double term2 = vt*std::erf(vt / sqrt4Dt) - (2 * sigma + vt)*std::erf((2 * sigma + vt) / sqrt4Dt);
   return 1.0 / 2 * (sqrt4Dt / M_SQRTPI*term1 + term2 + 2 * sigma);
}

// --------------------------------------------------------------------------------------------------------------------------------

double I_bd_r_3D(double r, double sigma, double t, double D)
{
   double Dt = D * t;
   double Dt2 = Dt + Dt;
   double Dt4 = Dt2 + Dt2;
   double sqrtDt = std::sqrt(Dt);
   double sqrtDt4 = std::sqrt(Dt4);
   double sigmasq = sigma * sigma;

   double sigmacb = sigmasq * sigma;
   double rcb = gsl_pow_3(r);

   double rsigma = r * sigma;
   double rps_sq = gsl_pow_2(r + sigma);
   double rms_sq = gsl_pow_2(r - sigma);

   double term1 = -2.0 * sqrtDt / M_SQRTPI;
   double term2 = std::exp(-sigmasq / Dt) * (sigmasq - Dt2);
   double term3 = -std::exp(-rps_sq / Dt4) * (rms_sq + rsigma - Dt2);
   double term4 = std::exp(-rms_sq / Dt4) * (rps_sq - rsigma - Dt2);
   double term5 = -sigmasq * 3.0 + Dt2;

   double term6 = (sigmacb - rcb) * std::erf((r - sigma) / sqrtDt4);
   double term7 = -(sigmacb + sigmacb) * std::erf(sigma / sqrtDt);
   double term8 = (sigmacb + rcb) * std::erf((r + sigma) / sqrtDt4);

   return (term1 * (term2 + term3 + term4 + term5) /* + sigmasq + rsigma + rsigma - Dt2)   //expm1 */ + term6 + term7 + term8) / 6.0;
}

// --------------------------------------------------------------------------------------------------------------------------------

double I_bd_r_1D(double r, double sigma, double t, double D, double v)
{
   if (D == 0) return 0;

   double Dt = D * t;
   double Dt2 = Dt + Dt;
   double Dt4 = Dt2 + Dt2;
   double sqrt4Dt = std::sqrt(Dt4);
   const double vt = v*t;

   double smrmvt_sq = gsl_pow_2(sigma - r - vt);
   double sprpvt_sq = gsl_pow_2(sigma + r + vt);
   double twospvt_sq = gsl_pow_2(2 * sigma + vt);

   double temp1 = -std::exp(-smrmvt_sq / Dt4) + std::exp(-sprpvt_sq / Dt4);
   double temp2 = std::exp(-vt*vt / Dt4) - std::exp(-twospvt_sq / Dt4);
   double term1 = sqrt4Dt / M_SQRTPI*(temp1 + temp2);

   double term2 = vt*std::erf(vt / sqrt4Dt) - (2 * sigma + vt)*std::erf((2 * sigma + vt) / sqrt4Dt);
   double term3 = (r - sigma + vt)*std::erf((sigma - r - vt) / sqrt4Dt);
   double term4 = (r + sigma + vt)*std::erf((r + sigma + vt) / sqrt4Dt);

   return 1.0 / 2 * (term1 + term2 + term3 + term4);
}

// --------------------------------------------------------------------------------------------------------------------------------

double drawR_gbd_3D(double rnd, double sigma, double t, double D)
{
   double I = I_bd_3D(sigma, t, D);

   double target = rnd * I;
   auto f = [sigma, t, D, target](double r) { return I_bd_r_3D(r, sigma, t, D) - target; };
   gsl_lambda<decltype(f)> F(f);

   double low = sigma;
   double high = sigma + 10.0 * std::sqrt(6.0 * D * t);

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-18, 1e-12, "drawR_gbd_3D");
}

// --------------------------------------------------------------------------------------------------------------------------------

double drawR_gbd_1D(double rnd, double sigma, double t, double D, double v)
{
   double I = I_bd_1D(sigma, t, D, v);

   double target = rnd * I;
   auto f = [sigma, t, D, target, v](double r) { return I_bd_r_1D(r, sigma, t, D, v) - target; };
   gsl_lambda<decltype(f)> F(f);

   double low = sigma;
   double high = sigma + 100.0 * std::sqrt(2.0 * D * t);

   root_fsolver_wrapper solver;
   solver.set(&F, low, high);
   return solver.findRoot(1e-18, 1e-12, "drawR_gbd_1D");
}

// --------------------------------------------------------------------------------------------------------------------------------
