#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include "ShellID.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct ShellCreateUtils
{
   // --------------------------------------------------------------------------------------------------------------------------------

   // The shell_distance_check is a MatrixSpace each_neighbor_xx processor that finds the shortest distance to the 'neighbors' from the given point.
   template<typename TMatrixSpace>
   class shell_distance_check
   {
      typedef typename TMatrixSpace::const_iterator         const_iterator;

   public:

      enum class Construct
      {
         SINGLE,
         OTHER,
         SHARE5050,        // share available space 50%-50% between creating a shell and next to another INIT-shell
         DEFUSION,         // share available space D1-D2 between creating a shell and next to another INIT-shell
      };

      shell_distance_check(ShellID sid, const Vector3& point, double distance, Construct construct) : ignore1_(sid), ignore2_(sid), point_(point), distance_(distance), construct_(construct) { }

      shell_distance_check(ShellID sid1, ShellID sid2, const Vector3& point, double distance, Construct construct) : ignore1_(sid1), ignore2_(sid2), point_(point), distance_(distance), construct_(construct) { }

      void operator()(const_iterator i, const Vector3& offset)
      {
         typename const_iterator::reference item(*i);
         if (item.first == ignore1_ || item.first == ignore2_) return;

         double distance = distance_;
         const Shell& shell = item.second;

         if (shell.code() == Shell::Code::INIT)    // init shell? leave some room for shell creation, so don't claim the full distance to an init shell.
         {
            ASSERT(shell.shape() == Shell::Shape::SPHERE);  // type is always sphere
            double radius = shell.get_sphere().radius();    // shell radius is particle radius
            switch (construct_)
            {
            case Construct::SINGLE: radius *= GfrdCfg.MULTI_SHELL_FACTOR; break;         // construct a Single, use fake shell for this neighbor with MSF
            case Construct::OTHER: radius *= GfrdCfg.SINGLE_SHELL_FACTOR; break;         // construct a Other, use fake shell for this neighbor with SSF
            case Construct::SHARE5050:
            {
               double length = (point_ - (shell.position() + offset)).length();          // distance between centers
               if (length / 2.0 < radius * GfrdCfg.MULTI_SHELL_FACTOR)                   // other particle should have at least radius*MSF
                  radius = radius * GfrdCfg.MULTI_SHELL_FACTOR;
               else
                  radius = length / 2.0;                                                 // use half the available space
            } break;

            //case Construct::DEFUSION: 
            //{
            //   double D2 = 0.5; // sim_.get_domain(shell.did())->particle()->D();         // what about Pairs or Multi's ? 
            //   double Dtot = D1_ + D2;
            //   if (D1_==0.0 || D2==0.0)
            //      radius *= GfrdCfg.MULTI_SHELL_FACTOR;        // one of two does not defuse, use MSF for scale
            //   else
            //      radius = (point_ - (shell.position() + offset)).length() * D2 / Dtot;         // share available space according to D1/D2 ratio
            //} break;

            default: ASSERT(0); break;
            }

            Sphere fake(shell.position() + offset, radius);
            distance = fake.distance(point_);
         }

         else if (shell.shape() == Shell::Shape::SPHERE)
            distance = shell.get_sphere().offset(offset).distance(point_);

         else if (shell.shape() == Shell::Shape::CYLINDER)
            distance = shell.get_cylinder().offset(offset).distance(point_);

         distance_ = std::min(distance_, distance / GfrdCfg.SAFETY);
      }

      double distance() const noexcept { return distance_; }

   private:
      const ShellID              ignore1_;   // ignore shells (its our own shell, when constructing for single, 
      const ShellID              ignore2_;   // or the shells from the two singles when forming a Pair.
      const Vector3&             point_;     // The point whose neighbors are being checked.
      double                     distance_;  // available space for shell creation
      const Construct            construct_;
   };

   // --------------------------------------------------------------------------------------------------------------------------------

   // The shell_interaction_check is a MatrixSpace each_neighbor_xx processor that finds particles (DomainID) that are within the interaction horizon
   template<typename TMatrixSpace>
   class shell_interaction_check
   {
      typedef typename TMatrixSpace::const_iterator         const_iterator;

   public:
      shell_interaction_check(ShellID sid, const Particle& particle) : ignore_(sid), point_(particle.position()), radius_(particle.radius()), hit_did_(0), multiple_() { }

      void operator()(const_iterator i, const Vector3& offset)
      {
         typename const_iterator::reference item(*i);
         if (item.first == ignore_) return;

         const Shell& shell = item.second;
         if (shell.code() != Shell::Code::INIT && shell.code() != Shell::Code::MULTI) return;      // look at init shells only (and multi, for new multi construction)

         ASSERT(shell.shape() == Shell::Shape::SPHERE);  // INIT shels and MULTI are always spherical
         double radius = shell.get_sphere().radius() / (shell.code() == Shell::Code::INIT ? 1.0 : GfrdCfg.MULTI_SHELL_FACTOR);     // InitShell size is particle, MultiShell size is MSF*particle!
         double distance = shell.shape() == Shell::Shape::SPHERE ? shell.get_sphere().offset(offset).distance(point_) : shell.get_cylinder().offset(offset).distance(point_);

         double pair_horizon = (radius + radius_) * GfrdCfg.SINGLE_SHELL_FACTOR;
         double multi_horizon = (radius + radius_) * GfrdCfg.MULTI_SHELL_FACTOR;

         // shell is init shell and particle is within single reaction horizon
         if (shell.code() == Shell::Code::INIT && distance < pair_horizon)
            hit_did_ = shell.did();

         // shell is init or multi shell and particle is within multi reaction horizon
         if (distance < multi_horizon)
         {
            auto did = shell.did();     // check if not already in list (can happen with multies, different shellID but same DomainID)
            if (multiple_.cend() == std::find(multiple_.cbegin(), multiple_.cend(), did))
               multiple_.emplace_back(shell.did());
         }
      }

      DomainID did() const noexcept { return hit_did_; }                     // when particle within pair_horizon
      bool multiple() const noexcept { return multiple_.size() > 1; }        // when two particles within multi_horizon
      gi::iteratorRange<std::vector<DomainID>> multi_range() const { return gi::iteratorRange<std::vector<DomainID>>(multiple_); }

   private:
      const ShellID              ignore_;   // structure (overlap checker) storing the overlapping particles.
      const Vector3&             point_;    // Particle location, the point whose neighbors are being checked.
      double                     radius_;   // Particle radius.

      DomainID                   hit_did_;   // when multiple is false, this contains zero (no hit) or the DomainID of (only) the particle within pair-interaction horizon (create a Pair domain)
      std::vector<DomainID>      multiple_;  // when multiple is true, there is more than one particle within the multi-interaction horizon
   };


   // --------------------------------------------------------------------------------------------------------------------------------

   // The shell_overlap_check is a MatrixSpace each_neighbor_xx processor that checks for overlapping shells (return domain list).
   template<typename TMatrixSpace>
   class shell_overlap_check
   {
      typedef typename TMatrixSpace::const_iterator         const_iterator;

   public:

      shell_overlap_check(const Sphere& sphere) : sphere_(sphere), overlap_domains_() { }

      void operator()(const_iterator i, const Vector3& offset)
      {
         const Shell& shell = (*i).second;

         if (shell.shape() == Shell::Shape::SPHERE)
         {
            double distance = shell.get_sphere().offset(offset).distance(sphere_.position());
            if (distance < sphere_.radius()) overlap_domains_.emplace_back(shell.did());
         }
         else if (shell.shape() == Shell::Shape::CYLINDER)
         {
            double distance = shell.get_cylinder().offset(offset).distance(sphere_.position());
            if (distance < sphere_.radius()) overlap_domains_.emplace_back(shell.did());
         }
      }

      std::vector<DomainID> overlap() { return std::move(overlap_domains_); }     // move, so only call once

   private:
      const Sphere&              sphere_;
      std::vector<DomainID>      overlap_domains_;
   };


   // --------------------------------------------------------------------------------------------------------------------------------

   // The shell_burst_check is a MatrixSpace each_neighbor_xx processor that shells within a given sphere, that may need bursting.
   template<typename TMatrixSpace>
   class shell_burst_check
   {
      typedef typename TMatrixSpace::const_iterator         const_iterator;

   public:

      shell_burst_check(const Sphere& sphere, std::vector<DomainID>& ignore) : ignore_(ignore), sphere_(sphere), burst_domains_() { }

      void operator()(const_iterator i, const Vector3& offset)
      {
         const Shell& shell = (*i).second;
         if (std::find(ignore_.cbegin(), ignore_.cend(), shell.did()) != ignore_.cend()) return;    // domain on ignore list?

         if (shell.code() == Shell::Code::INIT)
         {
            // If the domain was already bursted, but was not on the ignore list yet, put it there.
            // NOTE: we assume that the domain has already bursted around it when it was made!
            ignore_.emplace_back(shell.did());
         }
         else if (shell.shape() == Shell::Shape::SPHERE)
         {
            double distance = shell.get_sphere().offset(offset).distance(sphere_.position());
            if (distance < sphere_.radius()) burst_domains_.emplace_back(shell.did());
         }
         else if (shell.shape() == Shell::Shape::CYLINDER)
         {
            double distance = shell.get_cylinder().offset(offset).distance(sphere_.position());
            if (distance < sphere_.radius()) burst_domains_.emplace_back(shell.did());
         }
      }

      std::vector<DomainID> burst_range() { return std::move(burst_domains_); }     // move, so only call once

   private:
      std::vector<DomainID>&     ignore_;
      const Sphere&              sphere_;

      std::vector<DomainID>      burst_domains_;
   };

   // --------------------------------------------------------------------------------------------------------------------------------
};
