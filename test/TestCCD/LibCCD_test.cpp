// --------------------------------------------------------------------------------------------------------------------------------

#include <ccd/ccd.h>
#include "../../src/Common/ccdSupport.h"
#include "../../src/eGFRD/Vector3.hpp"
#include "LibCCD_test.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

int SphereSphereTest()
{
   ccd_t ccd;
   CCD_INIT(&ccd);
   ccd.support1 = ccdSupport;
   ccd.support2 = ccdSupport;
   //ccd.mpr_tolerance = 1e-9;
   //ccd.epa_tolerance = 1e-9;
   //ccd.dist_tolerance = 1e-9;

   const double scale = 1e-4;    // anything smaller does not work, regardless of tolerances, not sure why!

   ccd_sphere_t s1;
   s1.type = CCD_OBJ_SPHERE;
   s1.pos = { 8*scale, 8*scale, 0 };
   s1.quat = { 0., 0., 0., 1. };
   s1.radius = 5.9*scale;

   ccd_sphere_t s2;
   s2.type = CCD_OBJ_SPHERE;
   s2.pos = { 4*scale, 0, 0 };
   s2.quat = { 0., 0., 0., 1. };
   s2.radius = 1.15*scale;

   auto res1 = ccdGJKIntersect(&s1, &s2, &ccd);
   auto res2 = ccdGJKIntersect(&s2, &s1, &ccd);

   TINYTEST_EQUAL(res1, res2);
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int CylSphereCheck()
{
   ccd_t ccd;
   ccd_cyl_t c1, c2;

   c1.type = CCD_OBJ_CYL;
   c1.pos = { 0,0,0 };
   c1.quat = { 0,0,0,1 };

   ccd_vec3_t axis;
   ccdVec3Set(&axis, 0., 1., 0.);
   ccdQuatSetAngleAxis(&c1.quat, M_PI / 2., &axis);      // HORIZONTAL

   c2.type = CCD_OBJ_SPHERE;
   c2.pos = { 0,0,0 };
   c2.quat = { 0,0,0,1 };
   c2.radius = 1E-6;
   c2.height = 0;

   CCD_INIT(&ccd);
   ccd.support1 = ccdSupport;
   ccd.support2 = ccdSupport;

   c1.radius = 0.1;
   c1.height = 1.;      // TOTAL HEIGHT (halve on each side)

   ccdVec3Set(&c2.pos, 0.6, 0.0, 0.0);
   for (size_t i = 0; i < 60; i++)
   {
      int res = ccdGJKIntersect(&c1, &c2, &ccd);
      printf("%8.3lf\t%d\r\n", c2.pos.v[0], res);
      c2.pos.v[0] -= 0.01;
   }
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int OrientationCheck()
{
   Vector3 ori = Vector3::uz;
   Vector3 desired = Vector3::uz;
   auto axis = Vector3::cross(ori, desired);
   auto angle = std::acos(Vector3::dot(ori, desired) / (ori.length() * desired.length()));

   ccd_vec3_t a2 = ccd_vec3_t{ axis.X(), axis.Y(),axis.Z() };
   ccd_quat_t quat;
   ccdQuatSetAngleAxis(&quat, angle, &a2);


   ccd_t ccd;
   CCD_INIT(&ccd);
   ccd.support1 = ccdSupport;
   ccd.support2 = ccdSupport;

   ccd_sphere_t s;
   s.type = CCD_OBJ_SPHERE;
   s.pos = { 0, 0, 0 };
   s.quat = { 0., 0., 0., 1. };
   s.radius = 1e-03;

   ccd_cyl_t c;
   c.type = CCD_OBJ_CYL;
   c.pos = { 0, 0, 0 };
   c.quat = quat;
   c.radius = 0.2;
   c.height = 5;

   auto rIn = ccdGJKIntersect(&c, &s, &ccd);
   printf("%d=%s\r\n", rIn, "Inside");

   ccdVec3Set(&s.pos, -1, 0, 0);
   auto rLeft = ccdGJKIntersect(&c, &s, &ccd);
   printf("%d=%s\r\n", rLeft, "Left");
   ccdVec3Set(&s.pos, 1, 0, 0);
   auto rRight = ccdGJKIntersect(&c, &s, &ccd);
   printf("%d=%s\r\n", rRight, "Right");

   ccdVec3Set(&s.pos, 0, 0, -1);
   auto rBefore = ccdGJKIntersect(&c, &s, &ccd);
   printf("%d=%s\r\n", rBefore, "Before");
   ccdVec3Set(&s.pos, 0, 0, 1);
   auto rAfter = ccdGJKIntersect(&c, &s, &ccd);
   printf("%d=%s\r\n", rAfter, "After");

   ccdVec3Set(&s.pos, 0, -1, 0);
   auto rAbove = ccdGJKIntersect(&c, &s, &ccd);
   printf("%d=%s\r\n", rAbove, "Above");
   ccdVec3Set(&s.pos, 0, 1, 0);
   auto rBelow = ccdGJKIntersect(&c, &s, &ccd);
   printf("%d=%s\r\n", rBelow, "Below");
   
   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

int TestCylCyl(Vector3 ori, bool tX, bool tY, bool tZ)
{
   ccd_t ccd;
   CCD_INIT(&ccd);
   ccd.support1 = ccdSupport;
   ccd.support2 = ccdSupport;

   auto cyl = Cylinder(Vector3(), 0.2, ori, 2.5);
   auto sphere = Sphere(Vector3(), 1E-3);

   ccd_cyl_t c;
   ccdSetCyclinder(&c, cyl, 1);
   ccd_sphere_t s;
   ccdSetSphere(&s, sphere, 1);

   auto rIn = ccdGJKIntersect(&c, &s, &ccd);
   TINYTEST_ASSERT_MSG(rIn, "Inside");

   ccdVec3Set(&s.pos, -1, 0, 0);
   auto rLeft = ccdGJKIntersect(&c, &s, &ccd);
   TINYTEST_ASSERT_MSG(rLeft== (tX ? 1 : 0), "Left");
   ccdVec3Set(&s.pos, 1, 0, 0);
   auto rRight = ccdGJKIntersect(&c, &s, &ccd);
   TINYTEST_ASSERT_MSG(rRight== (tX ? 1 : 0), "Right");

   ccdVec3Set(&s.pos, 0, -1, 0);
   auto rAbove = ccdGJKIntersect(&c, &s, &ccd);
   TINYTEST_ASSERT_MSG(rAbove == (tY ? 1 : 0), "Above");
   ccdVec3Set(&s.pos, 0, 1, 0);
   auto rBelow = ccdGJKIntersect(&c, &s, &ccd);
   TINYTEST_ASSERT_MSG(rBelow == (tY ? 1 : 0), "Below");

   ccdVec3Set(&s.pos, 0, 0, -1);
   auto rBefore = ccdGJKIntersect(&c, &s, &ccd);
   TINYTEST_ASSERT_MSG(rBefore== (tZ ? 1 : 0), "Before");
   ccdVec3Set(&s.pos, 0, 0, 1);
   auto rAfter = ccdGJKIntersect(&c, &s, &ccd);
   TINYTEST_ASSERT_MSG(rAfter== (tZ ? 1 : 0), "After");
   
   return 1;
}

int OrientationTest()
{
   Vector3 ori[] = { Vector3::ux, Vector3::uy, Vector3::uz, };
   bool testX[] = { true, false, false, };
   bool testY[] = { false, true, false, };
   bool testZ[] = { false, false, true, };

   for (int i = 0; i < 3; ++i)
      TestCylCyl(ori[i], testX[i], testY[i], testZ[i]);

   return 1;
}

// --------------------------------------------------------------------------------------------------------------------------------

TINYTEST_START_SUITE(LibCCDTests);
TINYTEST_ADD_TEST(SphereSphereTest);
TINYTEST_ADD_TEST(CylSphereCheck);
TINYTEST_ADD_TEST(OrientationCheck);
TINYTEST_ADD_TEST(OrientationTest);
TINYTEST_END_SUITE();

// --------------------------------------------------------------------------------------------------------------------------------
