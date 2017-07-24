#include "Vector3.hpp"
#include "Matrix4.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT const Vector3 Vector3::null = Vector3();
GFRD_EXPORT const Vector3 Vector3::one = Vector3(1., 1., 1.);
GFRD_EXPORT const Vector3 Vector3::ux = Vector3(1., 0., 0.);
GFRD_EXPORT const Vector3 Vector3::uy = Vector3(0., 1., 0.);
GFRD_EXPORT const Vector3 Vector3::uz = Vector3(0., 0., 1.);


// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT Vector3 Vector3::transformVector(const Vector3& v, const Matrix4& mat)
{
   return mat.multiply(v);
}

// --------------------------------------------------------------------------------------------------------------------------------
