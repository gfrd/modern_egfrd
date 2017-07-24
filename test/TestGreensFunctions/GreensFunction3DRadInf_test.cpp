// --------------------------------------------------------------------------------------------------------------------------------

#include "GreensFunction3DRadInf_test.hpp"
#include "GreensFunction3DRadInf.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadInf_drawTime()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double r0 = 5e-8;
   auto gf = GreensFunction3DRadInf(D, kf, r0, sigma);

   double t = gf.drawTime(0.5);
   TINYTEST_FAIL(t <= 0.0);
   //TINYTEST_ALMOST_EQUALS(0.0008270075049176781, t);

   t = gf.drawTime(0.0);
   TINYTEST_FAIL(t < 0.0);
   //TINYTEST_ALMOST_EQUALS(2.666666666666667, t);

   t = gf.drawTime(0.9999999);
   TINYTEST_FAIL(t <= 0.0);
   //TINYTEST_ALMOST_EQUALS(1.1574074017561311e-05, t);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadInf_drawEventType()
{
   //double D = 1e-12;
   //double kf = 1e-8;
   //double sigma = 1e-8;
   //double a = 1e-7;
   //double r0 = 5e-8;
   //auto gf = GreensFunction3DRadInf(D, kf, r0, sigma, a);

   //double t = gf.drawTime(0.5);
   //auto event_type = gf.drawEventType(0.5, t);
   //TINYTEST_FAIL(event_type != GreensFunction::EventKind::IV_REACTION && event_type != GreensFunction::EventKind::IV_ESCAPE);

   //event_type = gf.drawEventType(0.0, t);
   //TINYTEST_EQUAL(GreensFunction::EventKind::IV_REACTION, event_type);

   //event_type = gf.drawEventType(0.999999, t);
   //TINYTEST_EQUAL(GreensFunction::EventKind::IV_ESCAPE, event_type);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadInf_drawR()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double r0 = 2e-8;

   auto gf = GreensFunction3DRadInf(D, kf, r0, sigma);

   double t = 1e-3;

   double r = gf.drawR(0.5, t);
   TINYTEST_FAIL(r < sigma);

   double r1 = gf.drawR(0.0, t);
   double r2 = gf.drawR(0.9999999, t);

   TINYTEST_FAIL(r1 < sigma);
   TINYTEST_FAIL(r2 < sigma);

   TINYTEST_FAIL(std::abs(r1 - sigma) > 1e-15);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadInf_drawTheta()
{
   double D = 1e-12;
   double kf = 1e-8;
   double sigma = 1e-8;
   double r0 = 2e-8;
   double t = 1e-3;
   double r = 2.1e-8;
   auto gf = GreensFunction3DRadInf(D, kf, r0, sigma);

   double theta = gf.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   theta = gf.drawTheta(0.0, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   theta = gf.drawTheta(0.9999999, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------


class UnitTestGreensFunction3DRadInf
{
public:
   static int test_p_int_r()
   {
      double D = 1e-12;
      double sigma = 1e-8;
      double kf = 1e-8;

      double t = 1e-3;
      double r0 = 2e-8;

      auto gf = GreensFunction3DRadInf(D, kf, r0, sigma);

      double pintr = gf.p_int_r(sigma, t);
      TINYTEST_EQUAL(0.0, pintr);

      return 1;
   }

   static int test_RnTable()
   {
      //double D = 1e-12;
      //double sigma = 1e-8;
      //double kf = 1e-8;
      //double a = 1e-7;
      //double r0 = 5e-8;

      //auto gf = GreensFunction3DRadInf(D, kf, r0, sigma, a);
      //double maxerror = 0.0;

      //for (int i = 0; i < 100; ++i)
      //{
      //   double alpha = gf.alpha0_i(i);
      //   double error = abs(gf.f_alpha0(alpha) / alpha);
      //   //#print error / alpha, gf.f_alpha0(alpha*1.1) / alpha
      //   maxerror = std::max(error, maxerror);
      //}

      //TINYTEST_FAIL(abs(maxerror) > 1e-10);
      return 1;
   }
};

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DRadInf_p_int_r()
{
   return UnitTestGreensFunction3DRadInf::test_p_int_r();
}
int GreensFunction3DRadInf_RnTable()
{
   return UnitTestGreensFunction3DRadInf::test_RnTable();
}

// --------------------------------------------------------------------------------------------------------------------------------


TINYTEST_START_SUITE(GreensFunction3DRadInf);
TINYTEST_ADD_TEST(GreensFunction3DRadInf_drawTime);
TINYTEST_ADD_TEST(GreensFunction3DRadInf_drawEventType)
TINYTEST_ADD_TEST(GreensFunction3DRadInf_drawR)
TINYTEST_ADD_TEST(GreensFunction3DRadInf_drawTheta)
TINYTEST_ADD_TEST(GreensFunction3DRadInf_p_int_r)
TINYTEST_ADD_TEST(GreensFunction3DRadInf_RnTable)
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
