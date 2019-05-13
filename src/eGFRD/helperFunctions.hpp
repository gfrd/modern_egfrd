#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <algorithm>
#include "Vector3.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

/**
* Return True if a and b are equal, subject to given tolerances.
*
* The (relative) tolerance must be positive and << 1.0
*
* Instead of specifying an absolute tolerance, you can specify a typical
* value for a or b. The absolute tolerance is then the relative tolerance
* multiplied by this typical value, and will be used when comparing a value
* to zero. By default, the typical value is 1.
*/
inline bool feq(double a, double b, double typical = 1., double tolerance = 1e-7)
{
    return std::abs(a - b) <= tolerance * (typical + std::min(std::abs(a), std::abs(b)));
}

// --------------------------------------------------------------------------------------------------------------------------------

/**
* Return True if a is greater than b, subject to given tolerances.
*
* The (relative) tolerance must be positive and << 1.0
*/
inline bool fgreater(double a, double b, double typical = 1., double tolerance = 1e-7)
{
    return a - b > tolerance * (typical + std::min(std::fabs(a), std::fabs(b)));
}

// --------------------------------------------------------------------------------------------------------------------------------


namespace cyclic
{

   /**
   * Transpose the position pos1 so that it can be used with another
   * position pos2.
   *
   * pos1 is transposed into one of mirror images of the cyclic boundary
   * condition so that the distance between pos1 and pos2 is smallest.
   *
   * Both of given pos1 and pos2 must be within the cyclic boundary.  However,
   * note that the returned transposed pos1 may not be within the cyclic boundary.
   */

   inline double cyclic_transpose(double p0, double p1, double world_size)
   {
      double diff(p1 - p0), half(world_size / 2);
      if (diff > half) return p0 + world_size;
      if (diff < -half)return p0 - world_size;
      return p0;
   }

   inline Vector3 cyclic_transpose(const Vector3& p0, const Vector3& p1, const Vector3& world_size)
   {
      return Vector3(cyclic_transpose(p0.X(), p1.X(), world_size.X()), cyclic_transpose(p0.Y(), p1.Y(), world_size.Y()), cyclic_transpose(p0.Z(), p1.Z(), world_size.Z()));
   }

   inline double distance_cyclic(const Vector3& p0, const Vector3& p1, const Vector3& world_size)
   {
      return (p0 - cyclic_transpose(p1, p0, world_size)).length();
   }

}

// --------------------------------------------------------------------------------------------------------------------------------
