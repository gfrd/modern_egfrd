// --------------------------------------------------------------------------------------------------------------------------------

#include <gsl/gsl_sf_bessel.h>
#include <algorithm>
#include "CylindricalBesselGenerator.hpp"
#include "CylindricalBesselGenerator_test.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// J: n range 0-50, size (11141..27853)
// Y: n range 0-50, size (22226..27797)

#if defined(_DEBUG)
#define Z_RANGE 500
#else
#define Z_RANGE 25000
#endif

// --------------------------------------------------------------------------------------------------------------------------------

int testCylindricalBessel0()
{
   const CylindricalBesselGenerator& cbt(CylindricalBesselGenerator::instance());
   const uint n = 0;

   for (uint i = 0; i < Z_RANGE; i++)
   {
      double z = 1 + i * 0.25;

      double j_table = cbt.J(n, z);
      double j_gsl = gsl_sf_bessel_J0(z);
      TINYTEST_ALMOST_EQUAL(j_table, j_gsl, 5.5E-8);

      double y_table = cbt.Y(n, z);
      double y_gsl = gsl_sf_bessel_Y0(z);
      TINYTEST_ALMOST_EQUAL(y_table, y_gsl, 5.5E-8);
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testCylindricalBessel1()
{
   const CylindricalBesselGenerator& cbt(CylindricalBesselGenerator::instance());
   const uint n = 1;

   for (uint i = 0; i < Z_RANGE; i++)
   {
      double z = 1 + i * 0.25;

      double j_table = cbt.J(n, z);
      double j_gsl = gsl_sf_bessel_J1(z);
      TINYTEST_ALMOST_EQUAL(j_table, j_gsl, 5.8E-8);

      double y_table = cbt.Y(n, z);
      double y_gsl = gsl_sf_bessel_Y1(z);
      TINYTEST_ALMOST_EQUAL(y_table, y_gsl, 5.5E-6);
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testCylindricalBesselN()
{
   const CylindricalBesselGenerator& cbt(CylindricalBesselGenerator::instance());

   for (uint n = cbt.MinNJ(); n < cbt.MaxNJ(); ++n)
   {
      for (uint i = 0; i < Z_RANGE; i++)
      {
         double z = 0.5 + i * 0.09;

         double j_table = cbt.J(n, z);
         double j_gsl = gsl_sf_bessel_Jn(n, z);

         double aerr = std::fabs(j_table - j_gsl);
         double rerr = aerr / j_gsl;

         TINYTEST_ASSERT(aerr < 5.8E-8 || rerr < 1E-12);       // abs.err can be large for huge values or rel.err can be large for small values (when either is error is small its OK).
      }
   }

   for (uint n = cbt.MinNY(); n < cbt.MaxNY(); ++n)
   {
      for (uint i = 0; i < Z_RANGE; i++)
      {
         double z = 0.5 + i * 0.09;

         double y_table = cbt.Y(n, z);
         double y_gsl = gsl_sf_bessel_Yn(n, z);

         double aerr = std::fabs(y_table - y_gsl);
         double rerr = aerr / y_gsl;

         TINYTEST_ASSERT(aerr < 3.1E-7 || rerr < 1E-12);
      }
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(CylindricalBesselGenerator);
TINYTEST_ADD_TEST(testCylindricalBessel0);
TINYTEST_ADD_TEST(testCylindricalBessel1);
TINYTEST_ADD_TEST(testCylindricalBesselN);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
