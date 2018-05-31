#include <sstream>
#include "Vector3.hpp"
#include "Vector3_test.hpp"
#include <randomNumberGenerator.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

class RandomNumberGenerator;

int testVector3Create()
{
   auto v1 = Vector3();
   TINYTEST_EQUAL_MSG(0.0, v1.X(), "default construct");
   TINYTEST_EQUAL_MSG(0.0, v1.Y(), "default construct");
   TINYTEST_EQUAL_MSG(0.0, v1.Z(), "default construct");

   auto v2 = Vector3(1, 2, 3);
   TINYTEST_EQUAL_MSG(1.0, v2.X(), "explicit construct");
   TINYTEST_EQUAL_MSG(2.0, v2.Y(), "explicit construct");
   TINYTEST_EQUAL_MSG(3.0, v2.Z(), "explicit construct");

   auto v3(v2);
   TINYTEST_EQUAL_MSG(1.0, v3.X(), "copy construct");
   TINYTEST_EQUAL_MSG(2.0, v3.Y(), "copy construct");
   TINYTEST_EQUAL_MSG(3.0, v3.Z(), "copy construct");

   v2 = v1;
   TINYTEST_EQUAL_MSG(0.0, v2.X(), "assignment operator");
   TINYTEST_EQUAL_MSG(0.0, v2.Y(), "assignment operator");
   TINYTEST_EQUAL_MSG(0.0, v2.Z(), "assignment operator");

   v1 = std::move(v3);
   TINYTEST_EQUAL_MSG(1.0, v1.X(), "assignment operator");
   TINYTEST_EQUAL_MSG(2.0, v1.Y(), "assignment operator");
   TINYTEST_EQUAL_MSG(3.0, v1.Z(), "assignment operator");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testVector3Normalize()
{
   auto v1 = Vector3(10, -20, 30);
   auto vn = v1.normal();

   TINYTEST_ALMOST_EQUAL_MSG(10 / 37.4165738677394138558374, vn.X(), 1E-15, "normalize");
   TINYTEST_ALMOST_EQUAL_MSG(-20 / 37.4165738677394138558374, vn.Y(), 1E-15, "normalize");
   TINYTEST_ALMOST_EQUAL_MSG(30 / 37.4165738677394138558374, vn.Z(), 1E-15, "normalize");

   TINYTEST_ALMOST_EQUAL_MSG(1.0, vn.length(), 1E-15, "normalize");

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testVector3Random()
{
   RandomNumberGenerator rng;
   auto v1 = Vector3::random(rng);
   TINYTEST_ALMOST_EQUAL_MSG(1.0, v1.length(), 1E-15, "random");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Vector3);
TINYTEST_ADD_TEST(testVector3Create);
TINYTEST_ADD_TEST(testVector3Normalize);
TINYTEST_ADD_TEST(testVector3Random);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
