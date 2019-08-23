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


namespace squared_distance
{
    /**
     * Calculates the minimal squared euclidian distance between two finite line segments in 3 dimensions.
     *
     * Say all points on segment 1, which begins at point start1 and ends at end1, can be described
     * as (start1 + lambda * u), where u is (end1 - start1), and lambda is a scalar between 0 and 1.
     * The same goes for segment 2: (start2 + mu * v).
     *
     * Now, we want to find those lambda and mu, for which the distance between the two points on the two
     * segments is minimal.
     *
     * Modified from http://geomalgorithms.com/a07-_distance.html#dist3D_Segment_to_Segment(), with
     * insights from Allen et al. (1993).
     * @param start1
     * @param end1
     * @param start2
     * @param end2
     * @return
     */
    inline double line_segments(Vector3 start1, Vector3 end1, Vector3 start2, Vector3 end2)
    {
        Vector3 u = end1 - start1;
        Vector3 v = end2 - start2;
        Vector3 w0 = start1 - start2;


        double a = Vector3::dot(u, u);
        double b = Vector3::dot(u, v);
        double c = Vector3::dot(v, v);
        double d = Vector3::dot(u, w0);
        double e = Vector3::dot(v, w0);
        double D = a*c - b*b;

        /*
         * (N is numerator, D is denominator)
         *
         * lambda = (be - cd) / (ac - b^2)
         *        = (be - cd) / D
         *        = N_lambda / D_lambda
         *
         * mu     = (ae - bd) / (ac - b^2)
         *        = (ae - bd) / D
         *        = N_mu / D_mu

         * The divisor terms are replaced by D_lambda and D_mu, respectively, to account for parallel segments in
         * the following code. The numerator terms are replaced by N_lambda, and N_mu, for the same reason.
         */
        double lambda, mu, N_lambda, N_mu, D_lambda, D_mu;

        if(feq(D, 0.0))
        {
            // Segments are (almost) parallel, so we set lambda = 0 and calculate the distance from mu = e/c;
            N_lambda = 0.0; // force using start1 of first segment
            D_lambda = 1.0; // to prevent possible division by 0.0 later
            N_mu = e;
            D_mu = c;
        }
        else
        {
            N_lambda = (b*e - c*d);
            N_mu = (a*e - b*d);
            if(N_lambda < 0.0)
            {
                // lambda < 0 -> the lambda = 0 edge is visible
                N_lambda = 0.0;
                N_mu = e;
                D_mu = c;
            }
            else if(N_lambda > D_lambda)
            {
                // lambda > 1 -> lambda = 1 edge is visible
                N_lambda = D_lambda;
                N_mu = e + b;
                D_mu = c;
            }
        }

        if(N_mu < 0.0)
        {
            // mu < 0 -> mu = 0 edge is visible
            N_mu = 0.0;

            // recompute lambda for this edge
            if (-d < 0.0)
            {
                N_lambda = 0.0;
            }
            else if(-d > a)
            {
                N_lambda = D_lambda;
            }
            else
            {
                N_lambda = -d;
                D_lambda = a;
            }
        }
        else if (N_mu > D_mu)
        {
            N_mu = D_mu;

            // recompute lambda for this edge
            if ((-d + b) < 0.0)
            {
                N_lambda = 0.0;
            }
            else if((-d + b) > a)
            {
                N_lambda = D_lambda;
            }
            else
            {
                N_lambda = -d + b;
                D_lambda = a;
            }
        }

        // Finally do the division to get lambda and mu
        lambda = (feq(fabs(N_lambda), 0, 0.0, 1e-20) ? 0.0 : N_lambda / D_lambda);
        mu = (feq(fabs(N_mu), 0, 0.0, 1e-20) ? 0.0 : N_mu / D_mu);

        // w is the vector between the two points on the two segments, that has the/a minimal length
        Vector3 w = w0 + (lambda * u) - (mu * v);

        // Return squared distance so prevent having to perform costly square root calculation.
        return Vector3::dot(w, w);
    }

    inline double line_segment_point(Vector3 start1, Vector3 end1, Vector3 point)
    {
        Vector3 u = end1 - start1;
        Vector3 v = point;
        Vector3 w0 = start1 - point;

        double lambda = -Vector3::dot(w0, u) / Vector3::dot(u, u);

        // If lambda < 0 or > 1, clamp to 0 respectively 1
        lambda = lambda < 0.0 ? 0.0 : (lambda > 1.0 ? 1.0 : lambda);

        Vector3 w = w0 + lambda * u;
        return Vector3::dot(w, w);
    }
}