#include <functional>
#include <sstream>
#include "Box.hpp"
#include "Box_test.hpp"
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box constructors. */
int testBoxCreate()
{
   Vector3 position(2, 3, 4);
   Vector3 half_extent(5, 6, 7);
   auto ux(Vector3::ux);
   auto uy(Vector3::uy);
   auto uz(Vector3::uz);

   Box b1;
   TINYTEST_EQUAL_MSG(Vector3::null, b1.position(), "default construct");

   Box b2(position);
   TINYTEST_EQUAL(position, b2.position());

   Box b3(position, half_extent);
   TINYTEST_EQUAL(position, b3.position());
   TINYTEST_EQUAL(half_extent, b3.half_extent());

   Box b4(position, half_extent, uz, uy, ux);
   TINYTEST_EQUAL(position, b4.position());
   TINYTEST_EQUAL(half_extent, b4.half_extent());
   TINYTEST_EQUAL(uz, b4.unit_x());
   TINYTEST_EQUAL(uy, b4.unit_y());
   TINYTEST_EQUAL(ux, b4.unit_z());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the degrees of freedom for particle movement on a Box. */
int testBoxDof()
{
    Box b1;
    TINYTEST_EQUAL(3, b1.dof());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Test the Box equal and not equal operator overloading. */
int testBoxEquality()
{
    Vector3 position(0.25, -1 / 8, -9.25);
    Vector3 half_extent(-0.5, 2.3, 7);
    auto uz = Vector3::uz;

    Box b1;
    Box b2;
    Box b3(position);
    Box b4(position, half_extent);
    Box b5(position, half_extent, uz, uz, uz);

    TINYTEST_ASSERT(b1 == b2);

    TINYTEST_ASSERT(b1 != b3);
    TINYTEST_ASSERT(b2 != b3);
    TINYTEST_ASSERT(b3 != b4);
    TINYTEST_ASSERT(b4 != b5);

    TINYTEST_ASSERT(b1 != b4);
    TINYTEST_ASSERT(b1 != b5);
    return 1;
}


// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box normal component for the ToInternal function. */
int testBoxToInternal()
{
    Vector3 position(1, 1, 1);
    Vector3 half_extent(-0.5, 2.3, 7);

    Box b1(position, half_extent);
    Vector3 v1(3, 2, 1);
    Vector3 vn = b1.to_internal(v1);

    TINYTEST_EQUAL(Vector3(2, 1, 0), vn);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box projection point function. */
int testBoxProject_point()
{
    Vector3 position(1, 2, 1);
    Vector3 half_extent(-0.5, 2.3, 7);

    Box b1(position, half_extent);
    Vector3 pos(1, 1, 0);
    auto pp = b1.project_point(pos);

    TINYTEST_EQUAL(Vector3(0, 0, 0), pp.first);
    TINYTEST_EQUAL(std::make_pair(0.0, 1.0), pp.second);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box projection of a point on surface function. */
int testBoxProject_point_on_surface()
{
    Box b1;
    Vector3 position(0, 0, 1);
    auto ppos = b1.project_point_on_surface(position);

    TINYTEST_EQUAL(std::make_pair(Vector3(), std::make_pair(double(),double())), ppos);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box distance function. */
int testBoxDistance()
{
    Box b1;
    Vector3 v1(1, 0, 0);
    auto vn = b1.distance(v1);

    TINYTEST_ALMOST_EQUAL_MSG(0.5, vn, 1E-15, "distance");

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box deflect function. */
int testBoxDeflect()
{
    Box b1;
    Vector3 r0(2, 3, 4);
    Vector3 d(1, 2, 3);

    auto df = b1.deflect(r0, d);
    TINYTEST_EQUAL(df.first, r0 + d);
    TINYTEST_EQUAL(df.second, false);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box random_position function. */
int testBoxRandom_position()
{
    Box b1;
    RandomNumberGenerator rng;
    auto vn = b1.random_position(rng);
    TINYTEST_ASSERT(vn.X() >= b1.position().X() - b1.half_extent().X());
    TINYTEST_ASSERT(vn.X() <= b1.position().X() + b1.half_extent().X());
    TINYTEST_ASSERT(vn.Y() >= b1.position().Y() - b1.half_extent().Y());
    TINYTEST_ASSERT(vn.Y() <= b1.position().Y() + b1.half_extent().Y());
    TINYTEST_ASSERT(vn.Z() >= b1.position().Z() - b1.half_extent().Z());
    TINYTEST_ASSERT(vn.Z() <= b1.position().Z() + b1.half_extent().Z());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Box print to stream function. */
int testBoxPrint()
{
   Box b1;

   std::stringstream sstream;
   sstream << b1;
   TINYTEST_STR_EQUAL("Box{P=(0, 0, 0), HE=(0.5, 0.5, 0.5), U=x(1, 0, 0),y(0, 1, 0),z(0, 0, 1)}", sstream.str().c_str());

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Box);
TINYTEST_ADD_TEST(testBoxCreate);
TINYTEST_ADD_TEST(testBoxDof);
TINYTEST_ADD_TEST(testBoxEquality);
TINYTEST_ADD_TEST(testBoxToInternal);
TINYTEST_ADD_TEST(testBoxProject_point);
TINYTEST_ADD_TEST(testBoxProject_point_on_surface);
TINYTEST_ADD_TEST(testBoxDistance);
TINYTEST_ADD_TEST(testBoxDeflect);
TINYTEST_ADD_TEST(testBoxRandom_position);
TINYTEST_ADD_TEST(testBoxPrint);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
