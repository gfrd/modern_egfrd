// --------------------------------------------------------------------------------------------------------------------------------

#define _USE_MATH_DEFINES
#include <math.h>
#include "GreensFunction3D_test.hpp"
#include "GreensFunction3D.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_instantiation()
{
   double D = 1e-12;
   double r0 = 5e-8;
   auto gf = GreensFunction3D(D, r0);

   TINYTEST_EQUAL(gf.getD(), D);
   TINYTEST_EQUAL(gf.getr0(), r0);
   TINYTEST_STR_EQUAL(gf.type_name(), "GreensFunction3D");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_drawR()
{
   double D = 1e-12;
   double r0 = 5e-8;
   double t = 1e-3;
   auto gf = GreensFunction3D(D, r0);

   auto r = gf.drawR(0.5, t);
   TINYTEST_ASSERT(r >= 0.0);

   r = gf.drawR(0.0, t);
   TINYTEST_ASSERT(r >= 0.0);

   r = gf.drawR(1.0, t);
   TINYTEST_ASSERT(r >= 0.0);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_drawR_zerot_is_r0()
{
   double D = 1e-12;
   double r0 = 2e-8;
   double t = 0.0;
   auto gf = GreensFunction3D(D, r0);

   auto r = gf.drawR(0.5, t);
   TINYTEST_EQUAL(r, r0);

   r = gf.drawR(0.0, t);
   TINYTEST_EQUAL(r, r0);

   r = gf.drawR(1.0, t);
   TINYTEST_EQUAL(r, r0);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_drawR_smallt()
{
   double D = 1e-12;
   double r0 = 2e-8;
   double t = 1e-4;
   auto gf = GreensFunction3D(D, r0);

   while (t > 1e-60)
   {
      double r = gf.drawR(0.5, t);

      TINYTEST_FAIL(r < 0.0);
      r = gf.drawR(0.0, t);
      TINYTEST_FAIL(r < 0.0);
      r = gf.drawR(1.0, t);
      TINYTEST_FAIL(r < 0.0);
      r = gf.drawR(1e-2, t);
      TINYTEST_FAIL(r < 0.0);

      t *= 1e-3;
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_draw_theta()
{
   double D = 1e-12;
   double r0 = 5e-8;
   double t = 1e-4;
   auto gf = GreensFunction3D(D, r0);

   //auto r = gf.drawR(0.5, r0, t);
   double r = r0;

   double theta = gf.drawTheta(0.5, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   theta = gf.drawTheta(0.0, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   theta = gf.drawTheta(1.0, r, t);
   TINYTEST_FAIL(theta < 0.0 || theta > M_PI);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

class UnitTestGreenFunction3D
{
public:

   static int ip_r_infinity_is_one()
   {
      double D = 1e-12;
      double t = 1e-5;
      double r0 = 5e-8;
      auto gf = GreensFunction3D(D, r0);
      auto ip = gf.ip_r(INFINITY, t);
      TINYTEST_EQUAL(1.0, ip);
      return 1;
   }

   static int int_p_r_is_ip_r()
   {
      double D = 1e-12;
      double t = 1e-5;
      double r0 = 5e-8;
      double  r = 2.5e-8;
      auto gf = GreensFunction3D(D, r0);

      auto ip = gf.ip_r(0.0, t);
      TINYTEST_EQUAL(0.0, ip);

      double maxr = 1e-6;

      int resolution = 20;
      for (int i = 1; i < resolution; i++)
      {
         r = i * maxr / resolution;
         ip = gf.ip_r(r, t);

         //auto result = scipy.integrate.quad(gf.p_r, 0.0, r, args = (t, ));
         //np = result[0]
         //self.assertAlmostEqual(0.0, (np - ip) / ip)
      }
      return 1;
   }
};

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_ip_r_infinity_is_one()
{
   return UnitTestGreenFunction3D::ip_r_infinity_is_one();
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_int_p_r_is_ip_r()
{
   return UnitTestGreenFunction3D::int_p_r_is_ip_r();
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_int_p_theta_is_p_r()
{
   double   D = 1e-12;
   //double t = 1e-5;
   double r0 = 5e-8;
   //double r = 2.5e-8;
   auto gf = GreensFunction3D(D, r0);

   //   auto ip = gf.ip_theta(M_PI, r, t);
      //result = scipy.integrate.quad(gf.p_theta, 0.0, M_PI, args = (r, t))
      //np = result[0]

     // auto pr = gf.p_r(r, t) / (2 * M_PI * r * r);

      //TINYTEST_EQUAL(0.0, (pr - ip) / pr);
      //TINYTEST_EQUAL(0.0, (pr - np) / pr);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3D_int_p_theta_is_ip_theta()
{
   double D = 1e-12;
   //double t = 1e-3;
   double r0 = 5e-8;
   //double r = 2.5e-8;
   auto gf = GreensFunction3D(D, r0);

   //   auto ip = gf.ip_theta(0.0, r, t);
   //   TINYTEST_EQUAL(0.0, ip);

   int resolution = 20;
   for (int i = 1; i < resolution; i++)
   {
      //auto theta = i * M_PI / resolution;
      //auto ip = gf.ip_theta(theta, r, t);
      //result = scipy.integrate.quad(gf.p_theta, 0.0, theta, args = (r, t))
      //np = result[0]
      //self.assertAlmostEqual(0.0, (np - ip) / ip)
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(GreensFunction3D);
TINYTEST_ADD_TEST(GreensFunction3D_instantiation);
TINYTEST_ADD_TEST(GreensFunction3D_drawR);
TINYTEST_ADD_TEST(GreensFunction3D_drawR_zerot_is_r0);
TINYTEST_ADD_TEST(GreensFunction3D_drawR_smallt);
TINYTEST_ADD_TEST(GreensFunction3D_draw_theta);
TINYTEST_ADD_TEST(GreensFunction3D_ip_r_infinity_is_one);
TINYTEST_ADD_TEST(GreensFunction3D_int_p_r_is_ip_r);
TINYTEST_ADD_TEST(GreensFunction3D_int_p_theta_is_p_r);
TINYTEST_ADD_TEST(GreensFunction3D_int_p_theta_is_ip_theta);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------

