#include "Matrix4.hpp"
#include "Vector3.hpp"
#include "DefsEgfrd.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT const Matrix4 Matrix4::null = Matrix4();

GFRD_EXPORT const Matrix4 Matrix4::identity = Matrix4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT Vector3 Matrix4::multiply(const Vector3& v) const
{
   return multiply(*this, v);
}

GFRD_EXPORT Matrix4 Matrix4::createRotationA(double theta, const Vector3& axis)
{
   // Returns a 3D (4x4) rotation matrix that rotates vectors theta radials around
   // the given axis.
   auto ux = axis.x_; auto uy = axis.y_; auto uz = axis.z_;       // from Wiki: Rotation_matrix
   auto cost = std::cos(theta);
   auto sint = std::sin(theta);

   return Matrix4(cost + ux*ux*(1.0 - cost), ux*uy*(1 - cost) - uz*sint, ux*uz*(1 - cost) + uy*sint, 0,
      uy*ux*(1 - cost) + uz*sint, cost + uy*uy*(1 - cost), uy*uz*(1 - cost) - ux*sint, 0,
      uz*ux*(1 - cost) - uy*sint, uz*uy*(1 - cost) + ux*sint, cost + uz*uz*(1 - cost), 0,
      0, 0, 0, 1);
}

GFRD_EXPORT Matrix4 Matrix4::createRotationAB(const Vector3& a, const Vector3& b)
{
   // Returns a 3D (4x4) rotation matrix that rotates vector a to vector b, around the axis
   // of a vector orthogonal to the plane spanned by (a, b).
   auto dot = Vector3::dot(a, b);

   if (dot > 0.9999999) return Matrix4::identity;
   if (dot < -0.999999) return Matrix4(-1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, 1);

   auto v = Vector3::cross(a, b);
   auto angle = std::acos(dot);
   return Matrix4::createRotationA(angle, v.normal());
}

GFRD_EXPORT Matrix4 Matrix4::createTranslate(const Vector3& t)
{
   return Matrix4(
      1, 0, 0, t.x_,
      0, 1, 0, t.y_,
      0, 0, 1, t.z_,
      0, 0, 0, 1);
}

GFRD_EXPORT Matrix4 Matrix4::createScale(const Vector3& t)
{
   return Matrix4(
      t.x_, 0, 0, 0,
      0, t.y_, 0, 0,
      0, 0, t.z_, 0,
      0, 0, 0, 1);
}

GFRD_EXPORT Matrix4 Matrix4::createCrossProduct(const Vector3& a)
{
   return Matrix4(
      0, -a.z_, a.y_, 0,
      a.z_, 0, -a.x_, 0,
      -a.y_, a.x_, 0, 0,
      0, 0, 0, 1);
}

GFRD_EXPORT Matrix4 Matrix4::createTensorProduct(const Vector3& a, const Vector3& b)
{
   return Matrix4(
      a.x_ * b.x_, a.x_ * b.y_, a.x_ * b.z_, 0,
      a.y_ * b.x_, a.y_ * b.y_, a.y_ * b.z_, 0,
      a.z_ * b.x_, a.z_ * b.y_, a.z_ * b.z_, 0,
      0, 0, 0, 1);
}

GFRD_EXPORT Vector3 Matrix4::multiply(const Matrix4& mat, const Vector3& v)
{
   return Vector3(
      mat.m00*v.x_ + mat.m10*v.y_ + mat.m20*v.z_ + mat.m30,
      mat.m01*v.x_ + mat.m11*v.y_ + mat.m21*v.z_ + mat.m31,
      mat.m02*v.x_ + mat.m12*v.y_ + mat.m22*v.z_ + mat.m32);
}

GFRD_EXPORT Matrix4 Matrix4::multiply(const Matrix4& a, const Matrix4& b)
{
   return Matrix4(
      a.m00 * b.m00 + a.m01 * b.m10 + a.m02 * b.m20 + a.m03 * b.m30,
      a.m00 * b.m01 + a.m01 * b.m11 + a.m02 * b.m21 + a.m03 * b.m31,
      a.m00 * b.m02 + a.m01 * b.m12 + a.m02 * b.m22 + a.m03 * b.m32,
      a.m00 * b.m03 + a.m01 * b.m13 + a.m02 * b.m23 + a.m03 * b.m33,
      
      a.m10 * b.m00 + a.m11 * b.m10 + a.m12 * b.m20 + a.m13 * b.m30,
      a.m10 * b.m01 + a.m11 * b.m11 + a.m12 * b.m21 + a.m13 * b.m31,
      a.m10 * b.m02 + a.m11 * b.m12 + a.m12 * b.m22 + a.m13 * b.m32,
      a.m10 * b.m03 + a.m11 * b.m13 + a.m12 * b.m23 + a.m13 * b.m33,
      
      a.m20 * b.m00 + a.m21 * b.m10 + a.m22 * b.m20 + a.m23 * b.m30,
      a.m20 * b.m01 + a.m21 * b.m11 + a.m22 * b.m21 + a.m23 * b.m31,
      a.m20 * b.m02 + a.m21 * b.m12 + a.m22 * b.m22 + a.m23 * b.m32,
      a.m20 * b.m03 + a.m21 * b.m13 + a.m22 * b.m23 + a.m23 * b.m33,
      
      a.m30 * b.m00 + a.m31 * b.m10 + a.m32 * b.m20 + a.m33 * b.m30,
      a.m30 * b.m01 + a.m31 * b.m11 + a.m32 * b.m21 + a.m33 * b.m31,
      a.m30 * b.m02 + a.m31 * b.m12 + a.m32 * b.m22 + a.m33 * b.m32,
      a.m30 * b.m03 + a.m31 * b.m13 + a.m32 * b.m23 + a.m33 * b.m33
      );
}

// --------------------------------------------------------------------------------------------------------------------------------
