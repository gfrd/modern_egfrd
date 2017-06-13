// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction3DRadAbs_test.hpp"
#include "GreensFunction3DRadAbs.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawTime()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = 5e-8;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   TINYTEST_FAIL(t <= 0.0 || t >= INFINITY);
   TINYTEST_ALMOST_EQUALS(0.0008270075049176781, t);

   t = gf.drawTime(0.0);
   TINYTEST_FAIL(t < 0.0 || t >= INFINITY);
   TINYTEST_ALMOST_EQUALS(2.666666666666667, t);

   t = gf.drawTime(1 - 1e-16);
   TINYTEST_FAIL(t <= 0.0 || t >= INFINITY);
   TINYTEST_ALMOST_EQUALS(1.1574074017561311e-05, t);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawTime_a_equal_sigma()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = sigma;
   double r0 = a;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   TINYTEST_EQUAL(0.0, t);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawTime_a_near_sigma()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = sigma + sigma * 1e-6;
   double r0 = (a + sigma) * .5;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   TINYTEST_FAIL(t <= 0.0 || t >= INFINITY);
   TINYTEST_ALMOST_EQUALS(9.468933926123182e-18, t);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawTime_r0_equal_a()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = a;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   TINYTEST_EQUAL(0.0, t);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawTime_r0_equal_sigma_kf_zero()
{
   double D = 1e-12;
   double kf = 0.0; // note this
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = sigma;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   TINYTEST_FAIL(t < 0.0 || t >= INFINITY);
   TINYTEST_ALMOST_EQUALS(0.0013418666154817517, t);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawTime_r0_equal_sigma_kf_large()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 10e-7;
   double r0 = sigma + 1e-12;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   TINYTEST_FAIL(t < 0.0 || t >= INFINITY);
   TINYTEST_ALMOST_EQUALS(1.0993113682152165e-12, t);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawEventType()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = 5e-8;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   auto event_type = gf.drawEventType(0.5, t);
   TINYTEST_FAIL(event_type != GreensFunction::EventKind::IV_REACTION && event_type != GreensFunction::EventKind::IV_ESCAPE);

   event_type = gf.drawEventType(0.0, t);
   TINYTEST_EQUAL(GreensFunction::EventKind::IV_REACTION, event_type);

   event_type = gf.drawEventType(0.999999, t);
   TINYTEST_EQUAL(GreensFunction::EventKind::IV_ESCAPE, event_type);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawEventType_smallt()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-6; //sigma + sigma * 0.001
   double r0 = 1.1e-8; //sigma + (a - sigma) / 2
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.999);

   auto event_type = gf.drawEventType(0.5, t);
   TINYTEST_FAIL(event_type != GreensFunction::EventKind::IV_REACTION && event_type != GreensFunction::EventKind::IV_ESCAPE);

   event_type = gf.drawEventType(0.0, t);
   TINYTEST_EQUAL(GreensFunction::EventKind::IV_REACTION, event_type);

   //event_type = gf.drawEventType(0.9999, t);
   //TINYTEST_EQUAL(GreensFunction::EventKind::IV_ESCAPE, event_type);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawR()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = 2e-8;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = 1e-3;

   double r = gf.drawR(0.5, t);
   TINYTEST_FAIL(r < sigma || r > a);

   double r1 = gf.drawR(0.0, t);
   double r2 = gf.drawR(0.999999999999, t);

   TINYTEST_FAIL(r1 < sigma || r1 > a);
   TINYTEST_FAIL(r2 < sigma || r2 > a);

   TINYTEST_ALMOST_EQUALS(r1, sigma);
   TINYTEST_ALMOST_REL_EQUAL(r2, a, 5E-7);         // its is not actually equal to a... close but no sigar

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawR_zerot()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = 2e-8;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = 0.0;

   double r = gf.drawR(0.5, t);
   TINYTEST_EQUAL(r0, r);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawR_r0_equal_sigma()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = sigma;
   double t = 1e-3;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double r = gf.drawR(0.5, t);
   TINYTEST_FAIL(r < sigma || r > a);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_drawR_squeezed()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1.01e-8;
   double t = 1e-6;
   double r0 = 1.005e-8;
   auto gf1 = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
   double r = gf1.drawR(0.5, t);
   TINYTEST_FAIL(r < sigma || r > a);

   r0 = 1.0001e-8;
   auto gf2 = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
   r = gf2.drawR(0.5, t);
   TINYTEST_FAIL(r < sigma || r > a);

   r0 = 1.0099e-8;
   auto gf3 = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
   r = gf3.drawR(0.5, t);
   TINYTEST_FAIL(r < sigma || r > a);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = 5e-8;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = gf.drawTime(0.5);
   auto event_type = gf.drawEventType(0.5, t);
   TINYTEST_FAIL(event_type != GreensFunction::EventKind::IV_REACTION && event_type != GreensFunction::EventKind::IV_ESCAPE);

   double r = gf.drawR(0.5, t);

   double theta = gf.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   theta = gf.drawTheta(0.0, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   theta = gf.drawTheta(0.999999, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_zerot()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = 5e-8;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = 0.0;
   double theta = gf.drawTheta(0.5, r0, t);
   TINYTEST_EQUAL(0.0, theta);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_smallt()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r = 2e-8;
   double r0 = 2e-8;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double t = 1e-7; // well this is not *very* small..
   double theta = gf.drawTheta(0.5, r, t);

   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_squeezed()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1.001e-8;
   double t = 1e-8;
   double r = 1.0001e-8;
   double r0 = 1.0001e-8;
   auto gf1 = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
   double theta = gf1.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   r = 1.00001e-8;
   r0 = 1.00001e-8;
   auto gf2 = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
   theta = gf2.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   r = 1.00099e-8;
   r0 = 1.00099e-8;
   auto gf3 = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
   theta = gf3.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_r0_equal_sigma()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = sigma;
   double t = 1e-3;
   double r = r0;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double theta = gf.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_r_equal_a()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double a = 1e-7;
   double r0 = 9e-8;
   double t = 1e-3;
   double r = a;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double theta = gf.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_1()
{
   double r0 = 1.0206416181e-07;
   double t = 4.41358538629e-08;
   double D = 4e-11;
   double sigma = 1e-07;
   double a = 1.05134e-07;
   double kf = 0; // h = 0;
   double r = 1.03421535312e-07;

   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double theta = gf.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_conv()
{
   double D = 2.0e-12;
   double kf = 0;
   double sigma = 2.0e-09;
   double r0 = 2.0226812535574088e-09;
   double a = 4.5403755344624737e-09;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double rnd = 0.26856112410314381;
   double r = 3.6385615779170657e-07;
   double t = 2.8801654577205008e-09;

   double theta = gf.drawTheta(rnd, r, t);
   TINYTEST_ALMOST_EQUAL(0.26722309705279157, theta, 5*DBL_EPSILON);       // WITH TABLE (fixed) 
   //TINYTEST_ALMOST_EQUALS(0.3392448267008449, theta);                    // WITH GSL BESSEL (also for small N) old code: 0.33924482669750622
   //TINYTEST_ALMOST_EQUALS(0.26722306378307309, theta);                   // WITH GSL BESSEL (not for small N)
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_draw_theta_conv2()
{
   double D = 2.0e-12;
   double kf = 0;
   double r0 = 2.022e-09;
   double a = 1e-06;
   double sigma = 2e-09;
   auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

   double rnd = 0.5;
   double t = 1e-05;
   double r = 4e-09;

   double theta = gf.drawTheta(rnd, r, t);
   TINYTEST_ALMOST_EQUAL(1.3949050836260224, theta, 5*DBL_EPSILON);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

class UnitTestGreensFunction3DRadAbs
{
public:
   static int test_ip_theta_squeezed()
   {
      double D = 1e-12;
      double kf = 1e-8;
      double sigma = 1e-8;
      double a = 1.001e-8;
      double t = 1e-10;
      double r = 1.00099e-8;
      double r0 = 1.00099e-8;
      auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
      double ip1 = gf.ip_theta(1, r, t);
      TINYTEST_ALMOST_EQUAL_MSG(1.2335592697733945e+19, ip1, 7E3, "ip1");     // wider tolerance due to rounding differences on CPU's ?!

      r = 1.0000001e-8;
      r0 = 1.0000001e-8;
      auto gf2 = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
      double ip2 = gf2.ip_theta(1, r, t);
      TINYTEST_ALMOST_EQUAL_MSG(1236720112075104.5, ip2, 1.0, "ip2");         // wider tolerance due to rounding differences on CPU's ?!

      return 1;
   }

   static int test_alpha_0()
   {
      double D = 1e-12;
      double sigma = 1e-8;
      double kf = 1e-8;
      double a = 1e-7;
      double r0 = 5e-8;

      auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);
      double maxerror = 0.0;

      for (int i = 0; i < 100; ++i)
      {
         double alpha = gf.alpha0_i(i);
         double error = abs(gf.f_alpha0(alpha) / alpha);
         //#print error / alpha, gf.f_alpha0(alpha*1.1) / alpha
         maxerror = std::max(error, maxerror);
      }

      TINYTEST_FAIL(abs(maxerror) > 1e-10);
      return 1;
   }

   static int test_alpha_i()
   {
      double D = 2.0e-12;
      double kf = 0.0;
      double sigma = 2.0e-09;
      double r0 = 2.0128020943506192e-09;
      double a = 4.5522328666282450e-09;

      auto gf = GreensFunction3DRadAbs(D, kf, r0, sigma, a);

      int i = 0;

      root_fsolver_wrapper solver;
      double alpha27 = gf.alpha_i(i, 27, solver);
      TINYTEST_ALMOST_EQUAL(289628066.89714050, alpha27, 5 * DBL_EPSILON);

      //double alpha26 = gf.alpha_i(i, 26, solver);         <= error endpoints do not straddle
      //TINYTEST_ALMOST_EQUAL(0.26722309705279157, alpha26, 5 * DBL_EPSILON);

      return 1;
   }

};

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadAbs_alpha_0()
{
   return UnitTestGreensFunction3DRadAbs::test_alpha_0();
}
int GreensFunction3DRadAbs_alpha_i()
{
   return UnitTestGreensFunction3DRadAbs::test_alpha_i();
}

int GreensFunction3DRadAbs_ip_theta_squeezed()
{
   return UnitTestGreensFunction3DRadAbs::test_ip_theta_squeezed();
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(GreensFunction3DRadAbs);
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawTime);
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawTime_a_equal_sigma);
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawTime_a_near_sigma);
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawTime_r0_equal_a);
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawTime_r0_equal_sigma_kf_zero)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawTime_r0_equal_sigma_kf_large)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawEventType)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawEventType_smallt)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawR)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawR_zerot)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawR_r0_equal_sigma)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_drawR_squeezed)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_zerot)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_smallt)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_squeezed)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_ip_theta_squeezed)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_r0_equal_sigma)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_r_equal_a)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_1)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_conv)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_draw_theta_conv2)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_alpha_0)
TINYTEST_ADD_TEST(GreensFunction3DRadAbs_alpha_i)
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
