// --------------------------------------------------------------------------------------------------------------------------------

#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_sum.h>
#include <gsl/gsl_machine.h>
#include <gsl/gsl_math.h>
#include "GslSumLevinU_test.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

#define N 500

int check_trunc(double * t, double expected, const char * desc);
int check_full(double * t, double expected, const char * desc);

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_zeta2()
{
   double t[N];
   int n;

   const double zeta_2 = M_PI * M_PI / 6.0;

   /* terms for zeta(2) */

   for (n = 0; n < N; n++)
   {
      double np1 = n + 1.0;
      t[n] = 1.0 / (np1 * np1);
   }

   check_trunc(t, zeta_2, "zeta(2)");
   check_full(t, zeta_2, "zeta(2)");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_exp10()
{
   double t[N];
   double x, y;
   int n;

   /* terms for exp(10.0) */
   x = 10.0;
   y = exp(x);

   t[0] = 1.0;
   for (n = 1; n < N; n++)
   {
      t[n] = t[n - 1] * (x / n);
   }

   check_trunc(t, y, "exp(10)");
   check_full(t, y, "exp(10)");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_expm10()
{
   double t[N];
   double x, y;
   int n;

   /* terms for exp(-10.0) */
   x = -10.0;
   y = exp(x);

   t[0] = 1.0;
   for (n = 1; n < N; n++)
   {
      t[n] = t[n - 1] * (x / n);
   }

   check_trunc(t, y, "exp(-10)");
   check_full(t, y, "exp(-10)");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_log12()
{
   double t[N];
   double x, y;
   int n;

   /* terms for -log(1-x) */
   x = 0.5;
   y = -log(1 - x);
   t[0] = x;
   for (n = 1; n < N; n++)
   {
      t[n] = t[n - 1] * (x * n) / (n + 1.0);
   }

   check_trunc(t, y, "-log(1/2)");
   check_full(t, y, "-log(1/2)");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_mlog2()
{
   double t[N];
   double x, y;
   int n;

   /* terms for -log(1-x) */
   x = -1.0;
   y = -log(1 - x);
   t[0] = x;
   for (n = 1; n < N; n++)
   {
      t[n] = t[n - 1] * (x * n) / (n + 1.0);
   }

   check_trunc(t, y, "-log(2)");
   check_full(t, y, "-log(2)");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_asymp()
{
   double t[N];
   int n;

   double result = 0.192594048773;

   /* terms for an alternating asymptotic series */

   t[0] = 3.0 / (M_PI * M_PI);

   for (n = 1; n < N; n++)
   {
      t[n] = -t[n - 1] * (4.0 * (n + 1.0) - 1.0) / (M_PI * M_PI);
   }

   check_trunc(t, result, "asymptotic series");
   check_full(t, result, "asymptotic series");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_eulergamma()
{
   double t[N];
   int n;

   /* Euler's gamma from GNU Calc (precision = 32) */

   double result = 0.5772156649015328606065120900824;

   /* terms for Euler's gamma */

   t[0] = 1.0;

   for (n = 1; n < N; n++)
   {
      t[n] = 1 / (n + 1.0) + log(n / (n + 1.0));
   }

   check_trunc(t, result, "Euler's constant");
   check_full(t, result, "Euler's constant");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_eta12()
{
   double t[N];
   int n;

   /* eta(1/2) = sum_{k=1}^{\infty} (-1)^(k+1) / sqrt(k)

   From Levin, Intern. J. Computer Math. B3:371--388, 1973.

   I=(1-sqrt(2))zeta(1/2)
   =(2/sqrt(pi))*integ(1/(exp(x^2)+1),x,0,inf) */

   double result = 0.6048986434216305;  /* approx */

                                        /* terms for eta(1/2) */

   for (n = 0; n < N; n++)
   {
      t[n] = (n % 2 ? -1 : 1) * 1.0 / sqrt(n + 1.0);
   }

   check_trunc(t, result, "eta(1/2)");
   check_full(t, result, "eta(1/2)");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GslSumLevinU_zero()
{
   double t[N];
   int n;

   double result = 1.23;

   for (n = 0; n < N; n++)
   {
      t[n] = (n == 0) ? 1.23 : 0.0;
   }

   check_trunc(t, result, "1.23 + 0 + 0 + 0...");
   check_full(t, result, "1.23 + 0 + 0 + 0...");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int check_trunc(double * t, double expected, const char * desc)
{
   double sum_accel, prec;
   gsl_sum_levin_utrunc_workspace * w = gsl_sum_levin_utrunc_alloc(N);
   gsl_sum_levin_utrunc_accel(t, N, w, &sum_accel, &prec);
   //gsl_test_rel(sum_accel, expected, 1e-8, "trunc result, %s", desc);
   TINYTEST_ALMOST_EQUAL_MSG(expected, sum_accel, 1e-8, desc);

   /* No need to check precision for truncated result since this is not a meaningful number */
   gsl_sum_levin_utrunc_free(w);
   return 1;
}

int check_full(double * t, double expected, const char * desc)
{
   double sum_accel, err_est, sd_actual, sd_est;
   gsl_sum_levin_u_workspace * w = gsl_sum_levin_u_alloc(N);
   gsl_sum_levin_u_accel(t, N, w, &sum_accel, &err_est);
   //gsl_test_rel(sum_accel, expected, 1e-8, "full result, %s", desc);
   TINYTEST_ALMOST_EQUAL_MSG(expected, sum_accel, 1e-8, desc);

   sd_est = -log10(err_est / fabs(sum_accel) + GSL_DBL_EPSILON);
   sd_actual = -log10(DBL_EPSILON + fabs((sum_accel - expected) / expected));

   /* Allow one digit of slop */
   //gsl_test(sd_est > sd_actual + 1.0, "full significant digits, %s (%g vs %g)", desc, sd_est, sd_actual);
   TINYTEST_FAIL_MSG(sd_est > sd_actual + 1.0, desc);
   gsl_sum_levin_u_free(w);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(GslSumLevinU);
TINYTEST_ADD_TEST(GslSumLevinU_zeta2);
TINYTEST_ADD_TEST(GslSumLevinU_exp10);
TINYTEST_ADD_TEST(GslSumLevinU_expm10)
TINYTEST_ADD_TEST(GslSumLevinU_log12)
TINYTEST_ADD_TEST(GslSumLevinU_mlog2)
TINYTEST_ADD_TEST(GslSumLevinU_asymp)
TINYTEST_ADD_TEST(GslSumLevinU_eulergamma)
TINYTEST_ADD_TEST(GslSumLevinU_eta12)
TINYTEST_ADD_TEST(GslSumLevinU_zero)
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
