// --------------------------------------------------------------------------------------------------------------------------------

#include <gsl/gsl_sf_bessel.h>
#include <algorithm>
#include "SphericalBesselGenerator.hpp"
#include "SphericalBesselGenerator_test.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// j: n range 4-51, size (11141..28978)
// y: n range 3-40, size (22276)

#if defined(_DEBUG)
#define Z_RANGE 500
#else
#define Z_RANGE 25000
#endif

// --------------------------------------------------------------------------------------------------------------------------------

int testSphericalBessel0()
{
   const SphericalBesselGenerator& sbt(SphericalBesselGenerator::instance());
   const uint n = 0;

   for (uint i = 0; i < Z_RANGE; i++)
   {
      double z = 1 + i * 0.25;

      double j_table = sbt.j(n, z);
      double j_gsl = gsl_sf_bessel_j0(z);
      TINYTEST_ALMOST_EQUAL(j_table, j_gsl, 1E-15);

      double y_table = sbt.y(n, z);
      double y_gsl = gsl_sf_bessel_y0(z);
      TINYTEST_ALMOST_EQUAL(y_table, y_gsl, 1E-15);
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testSphericalBessel1()
{
   const SphericalBesselGenerator& sbt(SphericalBesselGenerator::instance());
   const uint n = 1;

   for (uint i = 0; i < Z_RANGE; i++)
   {
      double z = 1 + i * 0.25;

      double j_table = sbt.j(n, z);
      double j_gsl = gsl_sf_bessel_j1(z);
      TINYTEST_ALMOST_EQUAL(j_table, j_gsl, 1E-15);

      double y_table = sbt.y(n, z);
      double y_gsl = gsl_sf_bessel_y1(z);
      TINYTEST_ALMOST_EQUAL(y_table, y_gsl, 1E-15);
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testSphericalBessel2()
{
   const SphericalBesselGenerator& sbt(SphericalBesselGenerator::instance());
   const uint n = 2;

   for (uint i = 0; i < Z_RANGE; i++)
   {
      double z = 1 + i * 0.25;

      double j_table = sbt.j(n, z);
      double j_gsl = gsl_sf_bessel_j2(z);
      TINYTEST_ALMOST_EQUAL(j_table, j_gsl, 1E-15);

      double y_table = sbt.y(n, z);
      double y_gsl = gsl_sf_bessel_y2(z);
      TINYTEST_ALMOST_EQUAL(y_table, y_gsl, 1E-15);
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testSphericalBessel3()
{
   const SphericalBesselGenerator& sbt(SphericalBesselGenerator::instance());
   const uint n = 3;

   for (uint i = 0; i < Z_RANGE; i++)
   {
      double z = 1 + i * 0.25;

      double j_table = sbt.j(n, z);
      double j_gsl = gsl_sf_bessel_jl(n, z);

      double y_table = sbt.y(n, z);
      double y_gsl = gsl_sf_bessel_yl(n, z);

      TINYTEST_ALMOST_EQUAL(j_table, j_gsl, 5E-15);
      TINYTEST_ALMOST_EQUAL(y_table, y_gsl, 5E-15);
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testSphericalBesselN()
{
   const SphericalBesselGenerator& sbt(SphericalBesselGenerator::instance());

   for (uint n = sbt.MinNj(); n < sbt.MaxNj(); ++n)
   {
      for (uint i = 0; i < Z_RANGE; i++)
      {
         double z = 0.5 + i * 0.09;

         double j_table = sbt.j(n, z);
         double j_gsl = gsl_sf_bessel_jl(n, z);

         double aerr = std::fabs(j_table - j_gsl);
         double rerr = aerr / j_gsl;

         TINYTEST_ASSERT(aerr < 5.7E-9 || rerr < 1E-12);       // abs.err can be large for huge values or rel.err can be large for small values (when either is error is small its OK).
         TINYTEST_ASSERT(aerr < 5.7E-9);
      }
   }

   for (uint n = sbt.MinNy(); n < sbt.MaxNy(); ++n)
   {
      for (uint i = 0; i < Z_RANGE; i++)
      {
         double z = 0.5 + i * 0.09;

         double y_table = sbt.y(n, z);
         double y_gsl = gsl_sf_bessel_yl(n, z);

         double aerr = std::fabs(y_table - y_gsl);
         double rerr = aerr / y_gsl;

         TINYTEST_ASSERT(aerr < 9.7E-8 || rerr < 1E-12);
         TINYTEST_ASSERT(aerr < 9.7E-8);
      }
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

//int testSphericalBesselNsmallZ()
//{
//   const SphericalBesselGenerator& sbt(SphericalBesselGenerator::instance());
//
//   uint n = 4;
//   {
//      for (uint i = 0; i < Z_RANGE; i++)
//      {
//         double z = i * 0.01;
//
//         double j_table = sbt.j(n, z);
//         double j_gsl = gsl_sf_bessel_jl(n, z);
//
//         double aerr = std::fabs(j_table - j_gsl);
//         double rerr = aerr / j_gsl;
//         
//         //TINYTEST_ASSERT(aerr < 5.7E-9 || rerr < 1E-12);
//         //TINYTEST_ASSERT(aerr < 5.7E-9);
//      }
//   }
//
//   {
//      for (uint i = 1; i < Z_RANGE; i++)
//      {
//         double z = i * 0.01;
//
//         double y_table = sbt.y(n, z);
//         double y_gsl = gsl_sf_bessel_yl(n, z);
//
//         double aerr = std::fabs(y_table - y_gsl);
//         double rerr = aerr / y_gsl;
//
//         //TINYTEST_ASSERT(aerr < 9.7E-8 || rerr < 1E-12);
//         //TINYTEST_ASSERT(aerr < 9.7E-8);
//      }
//   }
//   return 1;
//}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(SphericalBesselGenerator);
TINYTEST_ADD_TEST(testSphericalBessel0);
TINYTEST_ADD_TEST(testSphericalBessel1);
TINYTEST_ADD_TEST(testSphericalBessel2);
TINYTEST_ADD_TEST(testSphericalBessel3);
TINYTEST_ADD_TEST(testSphericalBesselN);
//TINYTEST_ADD_TEST(testSphericalBesselNsmallZ);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
