#include <functional>
#include <sstream>
#include "Disk.hpp"
#include "Disk_test.hpp"
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Disk's constructors. */
int testDiskCreate()
{
    Disk d1;
    TINYTEST_EQUAL_MSG(Vector3::null, d1.position(), "default construct");

    Vector3 position(1.25, 1 / 3, -1.25);
    double radius(8.75);
    auto unit_z = Vector3::ux;

    Disk d2(position, radius, unit_z);
    TINYTEST_EQUAL_MSG(position, d2.position(), "construct");
    TINYTEST_EQUAL(radius, d2.radius());
    TINYTEST_EQUAL(unit_z, d2.unit_z());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the degrees of freedom for particle movement on a Disk. */
int testDiskDof()
{
    Disk d1;
    TINYTEST_EQUAL(0, d1.dof());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Test the Disk's equal and not equal operator overloading. */
int testDiskEquality()
{
    Vector3 position(0.25, -1 / 8, -9.25);
    double radius(1.75);
    auto unit_z(Vector3::ux);

    Disk d1(position, radius, unit_z);
    Disk d2(position, radius, unit_z);
    Disk d3;

    TINYTEST_ASSERT(d1 == d2);
    TINYTEST_ASSERT(d1 != d3);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Disk's normal component for the ToInternal function.*/
int testDiskToInternal()
{
    Vector3 position(1, 1, 1);
    double radius(1);
    auto unit_z(Vector3::ux);

    Disk d1(position, radius, unit_z);
    Vector3 v1(3, 2, 1);
    auto vn = d1.to_internal(v1);

    TINYTEST_EQUAL(Vector2(1, 2), vn);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Disk's projection point function. */
int testDiskProject_point()
{
    Vector3 position(0, 0, 1);
    double radius(1);
    auto unit_z(Vector3::uz);

    Vector3 pos(1, 1, 0);
    Disk d1(position, radius, unit_z);
    auto pp = d1.project_point(pos);

    TINYTEST_EQUAL(Vector3(1,1,1), pp.first);
    TINYTEST_EQUAL(std::make_pair(-1.0, sqrt(2.0) - 1.0), pp.second);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Disk's projection of a point on surface function. */
int testDiskProject_point_on_surface()
{
    Vector3 position(0, 0, 1);
    double radius(1);
    auto unit_z = Vector3::uz;

    Vector3 pos(1, 1, 0);
    Disk d1(position, radius, unit_z);

    auto pp = d1.project_point(pos);
    auto ppos = d1.project_point_on_surface(pos);

    TINYTEST_EQUAL(pp, ppos);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Disks distance function. */
int testDiskDistance()
{
    Disk d1;
    Vector3 v1(1, 0, 0);
    auto vn = d1.distance(v1);

    TINYTEST_ALMOST_EQUAL_MSG(1.0, vn, 1E-15, "distance");
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Disks deflect function. */
int testDiskDeflect()
{
    Disk d1;
    Vector3 r0(2, 3, 4);
    Vector3 d(1, 2, 3);

    auto df = d1.deflect(r0, d);
    TINYTEST_EQUAL(df.first, r0 + d);
    TINYTEST_EQUAL(df.second, false);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Disk's random_position function. */
int testDiskRandom_position()
{
    Disk d1;
    RandomNumberGenerator rng;
    auto vn = d1.random_position(rng);

    TINYTEST_EQUAL(d1.position(), vn);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Disk);
TINYTEST_ADD_TEST(testDiskCreate);
TINYTEST_ADD_TEST(testDiskDof);
TINYTEST_ADD_TEST(testDiskEquality);
TINYTEST_ADD_TEST(testDiskToInternal);
TINYTEST_ADD_TEST(testDiskProject_point);
TINYTEST_ADD_TEST(testDiskProject_point_on_surface);
TINYTEST_ADD_TEST(testDiskDistance);
TINYTEST_ADD_TEST(testDiskDeflect);
TINYTEST_ADD_TEST(testDiskRandom_position);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
