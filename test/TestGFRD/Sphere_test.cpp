#include <functional>
#include <sstream>
#include "Sphere.hpp"
#include "Sphere_test.hpp"
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Sphere's constructors. */
int testSphereCreate()
{
    Sphere s1;
    TINYTEST_EQUAL_MSG(Vector3::null, s1.position(), "default construct");

    Vector3 position(1.25, 1/4, -1.25);
    double radius = 8.75;

    Sphere s2(position, radius);
    TINYTEST_EQUAL_MSG(position, s2.position(), "construct");
    TINYTEST_EQUAL(radius, s2.radius());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the degrees of freedom for particle movement on a Sphere. */
int testSphereDof()
{
    Sphere s1;
    TINYTEST_EQUAL(2, s1.dof());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Test the Sphere's equal and not equal operator overloading. */
int testSphereEquality()
{
    Vector3 position(0.25, -1 / 8, -9.25);
    double radius(1.75);

    Sphere s1(position, radius);
    Sphere s2(position, radius);
    Sphere s3;

    TINYTEST_ASSERT(s1 == s2);
    TINYTEST_ASSERT(s1 != s3);
    TINYTEST_ASSERT(s2 != s3);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Sphere's normal component for the ToInternal function.*/
int testSphereToInternal()
{
    Vector3 position(1, 1, 1);
    double radius(1);

    Sphere s1(position, radius);
    Vector3 v1(2, 2, 2);
    double si = s1.to_internal(v1);

    TINYTEST_EQUAL(sqrt(3.0), si);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Sphere's projection point function. */
int testSphereProject_point()
{
    Vector3 position(0, 0, 1);
    double radius(1);

    Sphere s1(position, radius);
    Vector3 pos(1, 1, 0);

    auto pp = s1.project_point(pos);

    TINYTEST_EQUAL(Vector3(0, 0, 1), pp.first);
    TINYTEST_EQUAL(std::make_pair(sqrt(3.0), 0.0), pp.second);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Sphere's projection of a point on surface function. */
int testSphereProject_point_on_surface()
{
    Vector3 position(0, 0, 1);
    double radius(1);

    Sphere s1(position, radius);
    Vector3 pos(1, 1, 0);
    auto ppos = s1.project_point_on_surface(pos);

    TINYTEST_EQUAL(std::make_pair(Vector3(), std::make_pair(double(), double())), ppos);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Spheres distance function. */
int testSphereDistance()
{
    Sphere s1;
    Vector3 v1(1, 0, 0);
    auto vn = s1.distance(v1);

    TINYTEST_ALMOST_EQUAL_MSG(1.0, vn, 1E-15, "distance");
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Spheres deflect function. */
int testSphereDeflect()
{
    Sphere s1;
    Vector3 r0(2, 3, 4);
    Vector3 d(1, 2, 3);

    auto df = s1.deflect(r0, d);
    TINYTEST_EQUAL(df.first, r0 + d);
    TINYTEST_EQUAL(df.second, false);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Sphere's random_position function. */
int testSphereRandom_position()
{
   Sphere s(Vector3(1,2,3), 5);

   RandomNumberGenerator rng;
   rng.seed(0xABCDEF);

   for (uint i = 0; i < 500; ++i)
   {
      auto vn = s.random_position(rng);           // point on the sphere surface
      double r = (vn - s.position()).length();
      TINYTEST_ASSERT_MSG(r <= s.radius(), "point not inside sphere!");
   }
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Sphere's print to stream function. */
int testSpherePrint()
{
    Sphere s1;

    std::stringstream sstream;
    sstream << s1;
    TINYTEST_STR_EQUAL("Sphere{P=(0, 0, 0), R=0}", sstream.str().c_str());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Sphere);
TINYTEST_ADD_TEST(testSphereCreate);
TINYTEST_ADD_TEST(testSphereDof);
TINYTEST_ADD_TEST(testSphereEquality);
TINYTEST_ADD_TEST(testSphereToInternal);
TINYTEST_ADD_TEST(testSphereProject_point);
TINYTEST_ADD_TEST(testSphereProject_point_on_surface);
TINYTEST_ADD_TEST(testSphereDistance);
TINYTEST_ADD_TEST(testSphereDeflect);
TINYTEST_ADD_TEST(testSphereRandom_position);
TINYTEST_ADD_TEST(testSpherePrint);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
