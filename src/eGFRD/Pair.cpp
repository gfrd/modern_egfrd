#include "Pair.hpp"
#include "exceptions.hpp"
#include "ShellCreateUtils.hpp"
#include "EGFRDSimulator.hpp"
#include <GreensFunction3DRadAbs.hpp>
#include <GreensFunction3DAbsSym.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

const double Pair::CUTOFF_FACTOR = 5.6;

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT bool PairSpherical::create_updated_shell(const shell_matrix_type& smat, const World& world, ShellID sid1, ShellID sid2)
{
   Vector3 pos1 = pid_pair1_.second.position();
   Vector3 pos2 = pid_pair2_.second.position();
   Vector3 pos2c = world.cyclic_transpose(pos2, pos1);

   Vector3 com = D_tot() > 0 ? (pid_pair2_.second.D() * pos1 + pid_pair1_.second.D() * pos2c) / D_tot() : 0.5 * (pos1 + pos2c);
   com_ = world.apply_boundary(com);

   iv_ = pos2c - pos1;

   THROW_UNLESS(no_space, r0() >= sigma());        // distance_from_sigma (pair gap) between %s and %s = %s < 0' % \(self.single1, self.single2, (self.r0 - self.sigma)))

   double min_radius = get_min_pair_size();
   double max_radius = smat.cell_size() * GfrdCfg.MAX_CELL_OCCUPANCY;      // assure spheres cannot overlap
   auto radius = max_radius;

   // when not diffusing, use minimal shell (we're not going anywhere, lease space for all others)
   if (D_tot() == 0.0) radius = min_radius;

   // check distances to surfaces, ignore def.struct and particle.structure
   for (auto s : world.get_structures())
   {
      if (/*s->id() != pid_pair_.second.structure_id() &&*/ s->id() != world.get_def_structure_id())    // ignore structure that particle is attached to
      {
         double distance = s->distance(com_);    // fix cyclic world !
         radius = std::min(radius, distance);
      }
   }

   shell_distance_checker sdc(sid1, sid2, com_, radius, shell_distance_checker::Construct::SHARE5050);
   CompileConfigSimulator::TBoundCondition::each_neighbor(smat, sdc, com_);
   radius = std::min(radius, sdc.distance());

   // OLDCODE radius = radius / 1.01;     // for compare with old code

   if (radius < min_radius) return false;             // no space for Pair domain.. it will be discarded, and a Multi is created instead!

   sid_pair_.second = Shell(domainID_, Sphere(com_, radius), Shell::Code::NORMAL);

   determine_radii(iv_.length(), radius);     // calculate inner and outer radius a_r_, a_R_.

   gf_com_ = std::make_unique<GreensFunction3DAbsSym>(D_R(), a_R_);
   gf_iv_ = std::make_unique<GreensFunction3DRadAbs>(D_tot(), k_total(), r0(), sigma(), a_r_);

   return true;
}

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT const PairGreensFunction& PairSpherical::choose_pair_greens_function() const
{
   double distance_from_sigma = r0() - sigma();
   double distance_from_shell = a_r_ - r0();
   double threshold_distance = CUTOFF_FACTOR * std::sqrt(6.0 * D_tot() * dt_);

   // if sigma reachable
   if (distance_from_sigma < threshold_distance)
   {
      // if shell reachable
      if (distance_from_shell < threshold_distance)

         // near both a and sigma;
         return *gf_iv_.get();

      // near sigma; use GreensFunction3DRadInf
      gf_tmp_ = std::make_unique<GreensFunction3DRadInf>(D_tot(), k_total(), r0(), sigma());
      return *gf_tmp_.get();
   }

   // sigma unreachable
   if (distance_from_shell < threshold_distance)
      // near a;
      gf_tmp_ = std::make_unique<GreensFunction3DAbs>(D_tot(), r0(), a_r_);
   else
      // distant from both a and sigma;
      gf_tmp_ = std::make_unique<GreensFunction3D>(D_tot(), r0());

   return *gf_tmp_.get();
}

// --------------------------------------------------------------------------------------------------------------------------------

GFRD_EXPORT void PairSpherical::determine_radii(double r0, double shell_size)
{
   // Determine a_r_ and a_R_ from the size of the protective domain.

   double radius1 = pid_pair1_.second.radius();
   double radius2 = pid_pair2_.second.radius();

   double D1 = pid_pair1_.second.D();
   double D2 = pid_pair2_.second.D();

   //double LD_MAX = INFINITY; // TODO was turned off (never greater than inf)
   //double a_r_max = LD_MAX * (r0 - sigma()) + sigma();
   // a_r_max is the maximum size of a_r_ for a given maximum ratio of l / delta

   // Make sure that D1 != 0 to avoid division by zero in the followings.
   //if D1 == 0:
//D1, D2 = D2, D1     // FIXME shouldn't we also here swap the particle radii if we're swapping diffusion contants ?
   //if __debug__ :
      //log.info('Swapping D1 and D2 because D1==0 (to avoid division by zero).')


   ASSERT(r0 >= sigma());        // '%s;  r0 %g < sigma %g' % (self, r0, self.sigma)

   double Da = D1;
   double Db = D2;
   double radiusa = radius1;
   double radiusb = radius2;

   // OLDCODE shell_size /= 1.01;     // for compare with old code

   // equalize expected mean t_r and t_R.
   if ((D_geom() - D2) * r0 / D_tot() + shell_size + std::sqrt(D2 / D1) * (radius1 - shell_size) - radius2 < 0)
   {
      Da = D2;
      Db = D1;
      radiusa = radius2;
      radiusb = radius1;
   }

   a_R_ = (D_geom() * (Db * (shell_size - radiusa) + Da * (shell_size - r0 - radiusa))) / (Da * Da + Da * Db + D_geom() * D_tot());
   a_r_ = (D_geom() * r0 + D_tot() * (shell_size - radiusa)) / (Da + D_geom());

   // Now if the planned domainsize for r is too large for proper convergence of the Greenfunctions, make it the maximum allowed size
   //if (a_r_ > a_r_max)
   //{
   //   a_r_ = a_r_max;
   //   a_R_ = shell_size - radiusa - a_r_ * Da / D_tot();
   //   Logger::get_logger("EGFRD").info("Domainsize changed for convergence: a_r_ = a_r_max = %.16g, a_R_ = %.16g", a_r_, a_R_);
   //}

   ASSERT(a_R_ + a_r_ * Da / D_tot() + radiusa >= a_R_ + a_r_ * Db / D_tot() + radiusb);
   ASSERT(std::abs(a_R_ + a_r_ * Da / D_tot() + radiusa - shell_size) < 1e-12 * shell_size);                        // here the shell_size is the relevant scale

   //   if __debug__:
   //log.debug('shell_size = %g, a_r_ = %g, a_R_ = %g r0 = %g' %
   //   (shell_size, a_r_, a_R_, r0))
   //   if __debug__ :
   //      tr = ((a_r_ - r0)**2) / (6 * self.D_r)       // the expected escape time of the IV
   //      if self.D_R == 0 :
   //         tR = numpy.inf
   //      else :
   //         tR = (a_R_**2) / (6 * self.D_R)          // the expected escape time of the CoM
   //         log.debug('tr = %g, tR = %g' % (tr, tR))


   ASSERT(a_r_ > 0);
   ASSERT(a_r_ > r0);
   ASSERT(a_R_ > 0 || (feq(a_R_, 0) && (D1 == 0 || D2 == 0)));
}


// --------------------------------------------------------------------------------------------------------------------------------
