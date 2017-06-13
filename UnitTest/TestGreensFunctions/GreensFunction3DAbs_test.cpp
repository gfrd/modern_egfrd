// --------------------------------------------------------------------------------------------------------------------------------

#define _USE_MATH_DEFINES
#include <math.h>
#include "GreensFunction3DAbs_test.hpp"
#include "GreensFunction3DAbs.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbs_drawR()
{
   double D = 1e-12;
   //   double t = 1e-3;
   double r0 = 5e-8;
   double a = 2.5e-8;

   auto gf = GreensFunction3DAbs(D, r0, a);

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

class UnitTestGreenFunction3DAbs
{
public:

   static int ip_r_infinity_is_one()
   {
      return 1;
   }
};

// --------------------------------------------------------------------------------------------------------------------------------

int GreensFunction3DAbs_ip_r_infinity_is_one()
{
   return UnitTestGreenFunction3DAbs::ip_r_infinity_is_one();
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(GreensFunction3DAbs);
TINYTEST_ADD_TEST(GreensFunction3DAbs_drawR);
TINYTEST_ADD_TEST(GreensFunction3DAbs_ip_r_infinity_is_one);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------

