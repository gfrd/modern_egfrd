#ifndef REVERSIBLE_HPP
#define REVERSIBLE_HPP

#include "Simulation3P.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Reversible : public Simulation3P
{
public:

   Reversible()
   {
      auto D = 1e-2;
      auto r = 2.5e-3;

      world_size_ = 1;
      Da_ = Db_ = Dc_ = D;
      ra_ = rb_ = rc_ = r;
      Na_ = Nb_ = 1; Nc_ = 0;
      kd_ = 100 * 2 * r * D;
      ds_ = 1e-9;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Reversible"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "ra" && args.isvalue(i + 1)) ra_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "Da" && args.isvalue(i + 1)) Da_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "rb" && args.isvalue(i + 1)) rb_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "Db" && args.isvalue(i + 1)) Db_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "ka" && args.isvalue(i + 1)) ka_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "kd" && args.isvalue(i + 1)) kd_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "ds" && args.isvalue(i + 1)) ds_ = args.value(++i);
      else return Simulation::HandleCommandArguments(i, args);
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------


   void print_usage() override
   {
      std::cout << "      -ra x.xx         Size of 'A' particle (radius in meter)\n";
      std::cout << "      -Da x.xx         Diffusion coefficient of particle 'A' in  m^2*s^-1\n";
      std::cout << "      -rb x.xx         Size of 'B' particle (radius in meter)\n";
      std::cout << "      -Db x.xx         Diffusion coefficient of particle 'B' in  m^2*s^-1\n";
      std::cout << "      -ka x.xx         Binding coefficient A+B particle in m^3*s^-1\n";
      std::cout << "      -kd x.xx         Unbinding coefficient C => A + B dissociation in  s^-1\n";
      std::cout << "      -ds x.xx         Separation between particle 'A' and 'B' (distance in meter)\n";
      Simulation::print_usage();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   bool SetupSimulation() override
   {
      sid_a = model_.add_species_type(SpeciesType("A", model_.get_def_structure_type_id(), Da_, ra_));
      sid_b = model_.add_species_type(SpeciesType("B", model_.get_def_structure_type_id(), Db_, rb_));
      sid_c = model_.add_species_type(SpeciesType("C", model_.get_def_structure_type_id(), Dc_, rc_));

      rules_.add_reaction_rule(ReactionRule(sid_c, kd_, std::vector<SpeciesTypeID>{sid_a, sid_b}));
      Simulation::SetupSimulation();

      auto pip_a = world_.add_particle(sid_a, world_.get_def_structure_id(), Vector3(0, 0, 0));
      pid_a_ = pip_a.first;

      pos_b_ = Vector3(rb_ + ra_ + ds_, 0, 0);
      auto pip_b = world_.add_particle(sid_b, world_.get_def_structure_id(), pos_b_);
      pid_b_ = pip_b.first;

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void SimStep(bool prerun) override
   {
      if (prerun) return; // ignore equalization phase

      Vector3 pos_b = world_.get_particle(pid_b_).second.position();
      total_dist_ += (pos_b_ - pos_b).length();

      if (simulator_->world().get_particle_ids(sid_c).size() > 0)
      {
         end_time_ = simulator_->time();     // stop!
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PostSimulation() override
   {
      std::cout << std::setw(14) << "time = " << simulator_->time() << " [s]\n";
      std::cout << std::setw(14) << "distance = " << total_dist_ << " [m]\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   double ds_;
   double total_dist_ = 0.0;        // traveled distance for particle B until it hit A and became a C particle
   Vector3 pos_b_;
   ParticleID pid_a_;
   ParticleID pid_b_;

};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /*REVERSIBLE_HPP*/