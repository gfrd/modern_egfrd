#ifndef GREENFUNCTIONHELPERS
#define GREENFUNCTIONHELPERS

// --------------------------------------------------------------------------------------------------------------------------------

//#define LOGALLGFDRAWS

#if defined(LOGALLGFDRAWS)
#define LOGDRAW(exp) exp;
#include "Logger.hpp"
#else
#define LOGDRAW(exp)
#endif

// --------------------------------------------------------------------------------------------------------------------------------

#include <GreensFunction.hpp>
#include <PairGreensFunction.hpp>
#include <cmath>
#include "randomNumberGenerator.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class GreenFunctionHelper
{
public:

   static double draw_time_wrapper(RandomNumberGenerator& rng, const GreensFunction& gf)
   {
      if (gf.getD() == 0) return INFINITY;
      double rnd = rng.uniform(0, 1);
      double t = gf.drawTime(rnd);
      LOGDRAW(Logger("CHECK").info("drawTime: %s, rnd=%.16g, %s => t=%.16g", gf.type_name(), rnd, gf.dump().c_str(), t));
      return t;
   }

   static GreensFunction::EventKind draw_eventtype_wrapper(RandomNumberGenerator& rng, const GreensFunction& gf, double dt)
   {
      double rnd = rng.uniform(0, 1);
      auto ek = gf.drawEventType(rnd, dt);
      LOGDRAW(Logger("CHECK").info("drawEventType: %s, rnd=%.16g, dt=%.16g, %s => ek=%s", gf.type_name(), rnd, dt, gf.dump().c_str(), ek==GreensFunction::EventKind::IV_ESCAPE?"ESCAPE":"REACTION"));
      return ek;
   }

   static double draw_reaction_time(double k, RandomNumberGenerator& rng)
   {
      if (k <= 0) return INFINITY;
      if (!std::isfinite(k)) return 0.0;
      double rnd = rng.uniform(0.0, 1.0);
      double t = rnd != 0 ? 1.0 / k * -std::log(rnd) : INFINITY;
      LOGDRAW(Logger("CHECK").info("drawReactionTime: rnd=%.16g, k=%.16g => t=%.16g", rnd, k, t));
      return t;
   }

   static double draw_r_wrapper(RandomNumberGenerator& rng, const GreensFunction& gf, double dt, double a, double sigma = -1)
   {
      double r;
      do
      {
         double rnd = rng.uniform(0, 1);
         r = gf.drawR(rnd, dt);
         LOGDRAW(Logger("CHECK").info("drawR: %s, rnd=%.16g, dt=%.16g, %s => r=%.16g", gf.type_name(), rnd, dt, gf.dump().c_str(), r));
      } while (r > a || (sigma >= 0 && r <= sigma));   // redraw; shouldn't happen often
      return r;
   }

   static double draw_theta_wrapper(RandomNumberGenerator& rng, const PairGreensFunction& gf, double r, double dt)
   {
      double rnd = rng.uniform(0, 1);
      double theta = gf.drawTheta(rnd, r, dt);              // theta [0,pi]
      LOGDRAW(Logger("CHECK").info("drawTheta: %s, rnd=%.16g, r=%.16g, dt=%.16g, %s => theta=%.16g", gf.type_name(), rnd, r, dt, gf.dump().c_str(), theta));
      return theta * (rng.uniform(0, 1) > 0.5 ? -1 : +1);     // between -pi and +pi
   }
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif

