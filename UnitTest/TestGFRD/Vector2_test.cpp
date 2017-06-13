#include <sstream>
#include "Vector2.hpp"
#include "Vector2_test.hpp"
#include <randomNumberGenerator.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

int testVector2Create()
{
    auto v1 = Vector2();
    TINYTEST_EQUAL_MSG(0.0, v1.X(), "default construct");
    TINYTEST_EQUAL_MSG(0.0, v1.Y(), "default construct");

    auto v2 = Vector2(1, 2);
    TINYTEST_EQUAL_MSG(1.0, v2.X(), "explicit construct");
    TINYTEST_EQUAL_MSG(2.0, v2.Y(), "explicit construct");

    auto v3 = v2;
    TINYTEST_EQUAL_MSG(1.0, v3.X(), "copy construct");
    TINYTEST_EQUAL_MSG(2.0, v3.Y(), "copy construct");

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testVector2Print()
{
    auto v1 = Vector2(10, -20);

    std::stringstream sstream;
    sstream << v1;
    TINYTEST_STR_EQUAL_MSG("(10, -20)", sstream.str().c_str(), "print");

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testVector2Normalize()
{
    auto v1 = Vector2(10, -20);
    auto vn = v1.normal();

    TINYTEST_ALMOST_EQUAL_MSG(1 / sqrt(5.0), vn.X(), 1E-15, "normalize");
    TINYTEST_ALMOST_EQUAL_MSG(-2 / sqrt(5.0), vn.Y(), 1E-15, "normalize");
    TINYTEST_ALMOST_EQUAL_MSG(1.0, vn.length(), 1E-15, "normalize");

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int testVector2Random()
{
   RandomNumberGenerator rng;
   auto v1 = Vector2::random(rng);
   TINYTEST_ALMOST_EQUAL_MSG(1.0, v1.length(), 1E-15, "random");
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Vector2);
TINYTEST_ADD_TEST(testVector2Create);
TINYTEST_ADD_TEST(testVector2Print);
TINYTEST_ADD_TEST(testVector2Normalize);
TINYTEST_ADD_TEST(testVector2Random);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
