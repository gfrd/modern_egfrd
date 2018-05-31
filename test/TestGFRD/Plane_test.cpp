#include <functional>
#include <sstream>
#include "Plane.hpp"
#include "Plane_test.hpp"
#include "randomNumberGenerator.hpp"


// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Plane constructors. */
int testPlaneCreate()
{
    Vector3 position(2, 3, 4);
    auto ux = Vector3::ux;
    auto uy = Vector3::uy;
    auto uz = Vector3::uz;
    double half_lx(2.25);
    double half_ly(1.25);
    bool is_one_sided(true);

    Plane p1;
    TINYTEST_EQUAL_MSG(Vector3::null, p1.position(), "default construct");

    Plane p2(position, ux, uy, half_lx, half_ly, is_one_sided);
    TINYTEST_EQUAL(position, p2.position());
    TINYTEST_EQUAL(ux, p2.unit_x());
    TINYTEST_EQUAL(uy, p2.unit_y());
    TINYTEST_EQUAL(uz, p2.unit_z());
    TINYTEST_EQUAL(Vector2(half_lx, half_ly), p2.half_extent());
    TINYTEST_EQUAL(is_one_sided, p2.is_one_sided());

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the degrees of freedom for particle movement on a Plane. */
int testPlaneDof()
{
    Plane p1;
    TINYTEST_EQUAL(2, p1.dof());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Test the Plane equal and not equal operator overloading. */
int testPlaneEquality()
{
    Vector3 position(2, 3, 4);
    auto ux = Vector3::ux;
    auto uy = Vector3::uy;
    double half_lx(2.25);
    double half_ly(1.25);
    bool is_one_sided(true);

    Plane p1;
    Plane p2;
    Plane p3(position, ux, uy, half_lx, half_ly, is_one_sided);

    TINYTEST_ASSERT(p1 == p2);
    TINYTEST_ASSERT(p1 != p3);
    TINYTEST_ASSERT(p2 != p3);
    return 1;
}


// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Plane normal component for the ToInternal function. */
int testPlaneToInternal()
{
    Vector3 position(2, 3, 4);
    auto ux = Vector3::ux;
    auto uy = Vector3::uy;
    double half_lx(2.25);
    double half_ly(1.25);
    bool is_one_sided(true);

    Vector3 v1(3, 2, 1);
    Plane p1(position, ux, uy, half_lx, half_ly, is_one_sided);
    auto vn = p1.to_internal(v1);

    TINYTEST_EQUAL(Vector3(1, -1, -3), vn);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Plane projection point function. */
int testPlaneProject_point()
{
    Vector3 position(2, 3, 4);
    auto ux = Vector3::ux;
    auto uy = Vector3::uy;
    double half_lx(2.25);
    double half_ly(1.25);
    bool is_one_sided(true);

    Plane p1(position, ux, uy, half_lx, half_ly, is_one_sided);
    Vector3 pos(1, 1, 0);
    auto pp = p1.project_point(pos);

    TINYTEST_EQUAL(Vector3(1, 1, 4), pp.first);
    TINYTEST_EQUAL(std::make_pair(-4.0, 1.0), pp.second);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Plane projection of a point on surface function. */
int testPlaneProject_point_on_surface()
{
    Plane p1;
    Vector3 position(0, 0, 1);
    auto ppos = p1.project_point_on_surface(position);

    TINYTEST_EQUAL(std::make_pair(Vector3(), std::make_pair(double(1.0), double(-0.5))), ppos);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Plane distance function. */
int testPlaneDistance()
{
    Plane p1;
    Vector3 v1(1, 0, 0);
    auto vn = p1.distance(v1);

    TINYTEST_ALMOST_EQUAL_MSG(0.5, vn, 1E-15, "distance");

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Plane deflect function. */
int testPlaneDeflect()
{
    Plane p1;
    Vector3 r0(2, 3, 4);
    Vector3 d(1, 2, 3);

    auto df = p1.deflect(r0, d);
    TINYTEST_EQUAL(df.first, r0 + d);
    TINYTEST_EQUAL(df.second, false);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Plane random_position function. */
int testPlaneRandom_position()
{
    Plane p;

    RandomNumberGenerator rng;
    auto vn = p.random_position(rng);
    TINYTEST_EQUAL(0.0, vn.Z());
    TINYTEST_ASSERT(std::abs(vn.X()) < p.half_extent().X());
    TINYTEST_ASSERT(std::abs(vn.Y()) < p.half_extent().Y());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Plane);
TINYTEST_ADD_TEST(testPlaneCreate);
TINYTEST_ADD_TEST(testPlaneDof);
TINYTEST_ADD_TEST(testPlaneEquality);
TINYTEST_ADD_TEST(testPlaneToInternal);
TINYTEST_ADD_TEST(testPlaneProject_point);
TINYTEST_ADD_TEST(testPlaneProject_point_on_surface);
TINYTEST_ADD_TEST(testPlaneDistance);
TINYTEST_ADD_TEST(testPlaneDeflect);
TINYTEST_ADD_TEST(testPlaneRandom_position);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
