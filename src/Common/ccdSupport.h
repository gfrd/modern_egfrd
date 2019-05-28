/***
 * libccd
 * ---------------------------------
 * Copyright (c)2010 Daniel Fiser <danfis@danfis.cz>
 *
 *
 *  This file is part of libccd.
 *
 *  Distributed under the OSI-approved BSD License (the "License");
 *  see accompanying file BDS-LICENSE for details or see
 *  <http://www.opensource.org/licenses/bsd-license.php>.
 *
 *  This software is distributed WITHOUT ANY WARRANTY; without even the
 *  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *  See the License for more information.
 */

 // --------------------------------------------------------------------------------------------------------------------------------

#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <ccd/quat.h>
#include "../eGFRD/Sphere.hpp"
#include "../eGFRD/Cylinder.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

// --------------------------------------------------------------------------------------------------------------------------------

#define CCD_OBJ_BOX 1
#define CCD_OBJ_SPHERE 2
#define CCD_OBJ_CYL 3

// --------------------------------------------------------------------------------------------------------------------------------

#define __CCD_OBJ__ \
    int type; \
    ccd_vec3_t pos; \
    ccd_quat_t quat;

   struct _ccd_obj_t {
      __CCD_OBJ__
   };
   typedef struct _ccd_obj_t ccd_obj_t;

   struct _ccd_box_t {
      __CCD_OBJ__
         ccd_real_t x, y, z; //!< Lengths of box's edges
   };
   typedef struct _ccd_box_t ccd_box_t;

   struct _ccd_sphere_t {
      __CCD_OBJ__
         ccd_real_t radius;
   };
   typedef struct _ccd_sphere_t ccd_sphere_t;

   struct _ccd_cyl_t {
      __CCD_OBJ__
         ccd_real_t radius;
      ccd_real_t height;
   };
   typedef struct _ccd_cyl_t ccd_cyl_t;

   // --------------------------------------------------------------------------------------------------------------------------------
      
   /**
    * Returns supporting vertex via v.
    * Supporting vertex is fathest vertex from object in direction dir.
    */
   static void ccdSupport(const void* _obj, const ccd_vec3_t * _dir, ccd_vec3_t * v)
   {
      // Support function is made according to Gino van den Bergen's paper
      //  A Fast and Robust CCD Implementation for Collision Detection of
      //  Convex Objects

      ccd_obj_t* obj = (ccd_obj_t*)_obj;
      ccd_vec3_t dir;
      ccd_quat_t qinv;

      ccdVec3Copy(&dir, _dir);
      ccdQuatInvert2(&qinv, &obj->quat);

      ccdQuatRotVec(&dir, &qinv);
      if (obj->type == CCD_OBJ_BOX)
      {
         ccd_box_t* box = (ccd_box_t*)obj;
         ccdVec3Set(v, ccdSign(ccdVec3X(&dir)) * box->x * CCD_REAL(0.5), ccdSign(ccdVec3Y(&dir)) * box->y * CCD_REAL(0.5), ccdSign(ccdVec3Z(&dir)) * box->z * CCD_REAL(0.5));
      }
      else if (obj->type == CCD_OBJ_SPHERE)
      {
         ccd_sphere_t* sphere = (ccd_sphere_t*)obj;
         ccd_real_t len;

         len = ccdVec3Len2(&dir);
         if (len - CCD_EPS > CCD_ZERO)
         {
            ccdVec3Copy(v, &dir);
            ccdVec3Scale(v, sphere->radius / CCD_SQRT(len));
         }
         else
         {
            ccdVec3Set(v, CCD_ZERO, CCD_ZERO, CCD_ZERO);
         }
      }
      else if (obj->type == CCD_OBJ_CYL)
      {
         ccd_cyl_t* cyl = (ccd_cyl_t*)obj;
         ccd_real_t zdist, rad;

         zdist = dir.v[0] * dir.v[0] + dir.v[1] * dir.v[1];
         zdist = CCD_SQRT(zdist);
         if (ccdIsZero(zdist))
         {
            ccdVec3Set(v, CCD_ZERO, CCD_ZERO, ccdSign(ccdVec3Z(&dir)) * cyl->height * CCD_REAL(0.5));
         }
         else
         {
            rad = cyl->radius / zdist;
            ccdVec3Set(v, rad * ccdVec3X(&dir), rad * ccdVec3Y(&dir), ccdSign(ccdVec3Z(&dir)) * cyl->height * CCD_REAL(0.5));
         }
      }

      // transform support vertex
      ccdQuatRotVec(v, &obj->quat);
      ccdVec3Add(v, &obj->pos);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   /**
    * Returns center of object.
    */
   static void ccdObjCenter(const void* obj, ccd_vec3_t* center);

   // --------------------------------------------------------------------------------------------------------------------------------

   // eGFRD Sphere to CCD sphere
   static void ccdSetSphere(ccd_sphere_t* s, const Sphere& sphere, double scale)
   {
      auto pos = scale * sphere.position();
      s->type = CCD_OBJ_SPHERE;
      s->pos = { pos.X(), pos.Y(), pos.Z() };
      s->quat = { 0., 0., 0., 1. };
      s->radius = scale * sphere.radius();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // eGFRD Cylinder to CCD cylinder
   static void ccdSetCyclinder(ccd_cyl_t* c, const Cylinder& cyl, double scale)
   {
      auto pos = scale * cyl.position();
      c->type = CCD_OBJ_CYL;
      c->pos = { pos.X(), pos.Y(), pos.Z() };
      c->radius = scale * cyl.radius();
      c->height = scale * 2 * cyl.half_length();

      // get quaternation for cylinder orientation
      Vector3 ori = cyl.unit_z().normal();
      Vector3 desired = Vector3::uz;
      auto axis = Vector3::cross(ori, desired);
      auto angle = std::acos(Vector3::dot(ori, desired));
      ccd_vec3_t a2 = ccd_vec3_t{ axis.X(), axis.Y(),axis.Z() };
      ccdQuatSetAngleAxis(&c->quat, angle, &a2);
   }

   // --------------------------------------------------------------------------------------------------------------------------------


#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

// --------------------------------------------------------------------------------------------------------------------------------
