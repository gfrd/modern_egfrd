#include <functional>
#include <sstream>
#include <helperFunctions.hpp>
#include "Cylinder.hpp"
#include "Cylinder_test.hpp"
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Cylinder's constructors. */
int testCylinderCreate()
{
    Cylinder c1;
    TINYTEST_EQUAL_MSG(Vector3::null, c1.position(), "default construct");

    Vector3 position(1.25, 1/3, -1.25);
    double radius(8.75);
    double half_length(2.25);
    auto unit_z(Vector3::uz);
    Cylinder c2(position, radius, unit_z, half_length);

    TINYTEST_EQUAL_MSG(position, c2.position(), "construct");
    TINYTEST_EQUAL(radius, c2.radius());
    TINYTEST_EQUAL(half_length, c2.half_length());
    TINYTEST_EQUAL(unit_z, c1.unit_z());

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the degrees of freedom for particle movement on a Cylinder. */
int testCylinderDof()
{
    Cylinder d1;
    TINYTEST_EQUAL(1, d1.dof());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Test the Cylinder's equal and not equal operator overloading. */
int testCylinderEquality()
{
    Vector3 position(0.25, -1/8, -9.25);
    double radius(1.75);
    double half_length(3.22);

    Cylinder d1(position, radius, Vector3::uz, half_length);
    Cylinder d2(position, radius, Vector3::uz, half_length);
    Cylinder d3;

    TINYTEST_ASSERT(d1 == d2);
    TINYTEST_ASSERT(d1 != d3);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Cylinder's normal component for the ToInternal function.*/
int testCylinderToInternal()
{
    Vector3 position(1, 1, 1);
    double radius(1);
    double half_length(2);

    Cylinder d1(position, radius, Vector3::ux, half_length);
    Vector3 v1(3, 2, 1);
    Vector2 vn = d1.to_internal(v1);

    TINYTEST_EQUAL(Vector2(1, 2), vn);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Cylinder's projection point function. */
int testCylinderProject_point()
{
    Vector3 position(0, 0, 1);
    double radius(1);
    double half_length(3.0);

    Cylinder d1(position, radius, Vector3::ux, half_length);
    Vector3 pos(1, 1, 0);
    auto pp = d1.project_point(pos);

    TINYTEST_EQUAL(Vector3(sqrt(2), 0, 1), pp.first);
    TINYTEST_EQUAL(std::make_pair(sqrt(2.0), - 2.0), pp.second);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Cylinder's projection of a point on surface function. */
int testCylinderProject_point_on_surface()
{
    Vector3 position(1, 0, 1);
    double radius(1);
    double half_length(1.0);

    Cylinder d1(position, radius, Vector3::ux, half_length);
    Vector3 pos(1, 0, 0);

    auto ppos = d1.project_point_on_surface(pos);
    auto res = std::make_pair(Vector3(1, 0, 0), std::make_pair(0.0,-1.0));

    TINYTEST_EQUAL(res, ppos);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Cylinders distance function. */
int testCylinderDistance()
{
    Cylinder d1;
    Vector3 v1(1, 0, 1);
    auto vn = d1.distance(v1);

    TINYTEST_ALMOST_EQUAL_MSG(0.5, vn, 1E-15, "distance");
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Cylinders deflect function. */
int testCylinderDeflect()
{
    Cylinder d1;
    Vector3 r0(2, 3, 4);
    Vector3 d(1, 2, 3);

    auto df = d1.deflect(r0, d);
    TINYTEST_EQUAL(df.first, r0 + d);
    TINYTEST_EQUAL(df.second, false);
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests the Cylinder's random_position function. */
int testCylinderRandom_position()
{
    Cylinder c;

    RandomNumberGenerator rng;

    auto vn = c.random_position(rng);
    TINYTEST_ASSERT(vn.X() >= c.position().X() - c.half_length());      // point along length of cylinder
    TINYTEST_ASSERT(vn.X() <= c.position().X() + c.half_length());
    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests distance between cylinders. */
int testCylinderDistanceCylinderParallel()
{
    Cylinder c1(Vector3(0.0, 0.0, 0.0), 1.0, Vector3(0.0, 1.0, 0.0), 0.5),
            c2(Vector3(3.0, 0.0, 0.0), 1.0, Vector3(0.0, 1.0, 0.0), 0.5);

    Vector3 start1 = c1.position() - c1.half_length() * c1.unit_z(),
            end1 = c1.position() + c1.half_length() * c1.unit_z(),
            start2 = c2.position() - c2.half_length() * c2.unit_z(),
            end2 = c2.position() + c2.half_length() * c2.unit_z();

    double line_squared_dist = squared_distance::line_segments(start1, end1, start2, end2);
    double sum_radii = c1.radius() + c2.radius();

    TINYTEST_ASSERT(line_squared_dist == 3.0*3.0);
    TINYTEST_ASSERT(std::sqrt(line_squared_dist) - sum_radii == 1.0);

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests distance between cylinders. */
int testCylinderDistanceCylinderParallelOffset()
{
    Cylinder c1(Vector3(0.0, 0.0, 0.0), 1.0, Vector3(0.0, 1.0, 0.0), 0.5),
            c2(Vector3(3.0, 2.0, 0.0), 1.0, Vector3(0.0, 1.0, 0.0), 0.5);

    Vector3 start1 = c1.position() - c1.half_length() * c1.unit_z(),
            end1 = c1.position() + c1.half_length() * c1.unit_z(),
            start2 = c2.position() - c2.half_length() * c2.unit_z(),
            end2 = c2.position() + c2.half_length() * c2.unit_z();

    double line_squared_dist = squared_distance::line_segments(start1, end1, start2, end2);

    TINYTEST_ASSERT(line_squared_dist == 3.0*3.0);

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests distance between cylinders. */
int testCylinderDistanceCylinderOrthogonal()
{
    Cylinder c1(Vector3(0.0, 0.0, 0.0), 1.0, Vector3(0.0, 1.0, 0.0), 0.5),
            c2(Vector3(0.0, 0.0, 3.0), 1.0, Vector3(0.0, 0.0, 1.0), 0.5);

    Vector3 start1 = c1.position() - c1.half_length() * c1.unit_z(),
            end1 = c1.position() + c1.half_length() * c1.unit_z(),
            start2 = c2.position() - c2.half_length() * c2.unit_z(),
            end2 = c2.position() + c2.half_length() * c2.unit_z();

    double line_squared_dist = squared_distance::line_segments(start1, end1, start2, end2);
    double dist = std::sqrt(line_squared_dist);

    TINYTEST_ASSERT(feq(line_squared_dist, 2.5 * 2.5));
    TINYTEST_ASSERT(feq(dist, 2.5));

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests distance between a cylinder and point. */
int testCylinderDistancePointOrthogonal()
{
    Cylinder c1(Vector3(0.0, 0.0, 0.0), 1.0, Vector3(0.0, 1.0, 0.0), 0.5);

    Vector3 start1 = c1.position() - c1.half_length() * c1.unit_z(),
            end1 = c1.position() + c1.half_length() * c1.unit_z(),
            point = Vector3(0, 0.0, 2.0);

    double line_squared_dist = squared_distance::line_segment_point(start1, end1, point);

    TINYTEST_ASSERT(feq(line_squared_dist, 2.0 * 2.0));

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

/* Tests distance between a cylinder and point. */
int testCylinderDistancePointSkewed()
{
    Cylinder c1(Vector3(0.0, 0.0, 0.0), 1.0, Vector3(0.0, 1.0, 0.0), 0.5);

    Vector3 start1 = c1.position() - c1.half_length() * c1.unit_z(),
            end1 = c1.position() + c1.half_length() * c1.unit_z(),
            point = Vector3(0, 2.5, 2.0);

    double line_squared_dist = squared_distance::line_segment_point(start1, end1, point);

    TINYTEST_ASSERT(feq(line_squared_dist, 2*2 + 2*2));

    return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(Cylinder);
TINYTEST_ADD_TEST(testCylinderCreate);
TINYTEST_ADD_TEST(testCylinderDof);
TINYTEST_ADD_TEST(testCylinderEquality);
TINYTEST_ADD_TEST(testCylinderToInternal);
TINYTEST_ADD_TEST(testCylinderProject_point);
TINYTEST_ADD_TEST(testCylinderProject_point_on_surface);
TINYTEST_ADD_TEST(testCylinderDistance);
TINYTEST_ADD_TEST(testCylinderDeflect);
TINYTEST_ADD_TEST(testCylinderRandom_position);
TINYTEST_ADD_TEST(testCylinderDistanceCylinderParallel);
TINYTEST_ADD_TEST(testCylinderDistanceCylinderParallelOffset);
TINYTEST_ADD_TEST(testCylinderDistanceCylinderOrthogonal);
TINYTEST_ADD_TEST(testCylinderDistancePointOrthogonal);
TINYTEST_ADD_TEST(testCylinderDistancePointSkewed);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
