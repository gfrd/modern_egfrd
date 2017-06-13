#ifndef SIMULATION3P_HPP
#define SIMULATION3P_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "Simulation.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

// Represents a EGFRD simulation with reversible reaction A + B <--> C.
class Simulation3P : public Simulation
{
public:
   // --------------------------------------------------------------------------------------------------------------------------------

   explicit Simulation3P() : Simulation() { }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "Na" && args.isvalue(i + 1)) Na_ = args.number(++i);
      else if (args.isparam(i) && args.option(i) == "ra" && args.isvalue(i + 1)) ra_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "Da" && args.isvalue(i + 1)) Da_ = args.value(++i);

      else if (args.isparam(i) && args.option(i) == "Nb" && args.isvalue(i + 1)) Nb_ = args.number(++i);
      else if (args.isparam(i) && args.option(i) == "rb" && args.isvalue(i + 1)) rb_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "Db" && args.isvalue(i + 1)) Db_ = args.value(++i);

      else if (args.isparam(i) && args.option(i) == "Nc" && args.isvalue(i + 1)) Nc_ = args.number(++i);
      else if (args.isparam(i) && args.option(i) == "rc" && args.isvalue(i + 1)) rc_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "Dc" && args.isvalue(i + 1)) Dc_ = args.value(++i);

      else if (args.isparam(i) && args.option(i) == "ka" && args.isvalue(i + 1)) ka_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "kd" && args.isvalue(i + 1)) kd_ = args.value(++i);
      else return Simulation::HandleCommandArguments(i, args);
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      std::cout << "      -Na n            Initial number of 'A' particles" << "\n";
      std::cout << "      -ra x.xx         Size of 'A' particle (radius in meter)" << "\n";
      std::cout << "      -Da x.xx         Diffusion coefficient of particle 'A' in  m^2*s^-1" << "\n";
      std::cout << "      -Nb n            Initial number of 'B' particles" << "\n";
      std::cout << "      -rb x.xx         Size of 'B' particle (radius in meter)" << "\n";
      std::cout << "      -Db x.xx         Diffusion coefficient of particle 'B' in  m^2*s^-1" << "\n";
      std::cout << "      -Nc n            Initial number of 'C' particles" << "\n";
      std::cout << "      -rc x.xx         Size of 'C' particle (radius in meter)" << "\n";
      std::cout << "      -Dc x.xx         Diffusion coefficient of particle 'C' in  m^2*s^-1" << "\n";
      std::cout << "      -ka x.xx         Binding coefficient A+B particle in m^3*s^-1" << "\n";
      std::cout << "      -kd x.xx         Unbinding coefficient C => A + B dissociation in  s^-1" << "\n" << "\n";
      Simulation::print_usage();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   void PrintSettings() override
   {
      Simulation::PrintSettings();

      std::cout << std::setw(14) << "#A = " << Na_ << "\n";
      std::cout << std::setw(14) << "D_a = " << Da_ << " [m^2*s^-1]" << "\n";
      std::cout << std::setw(14) << "r_a = " << ra_ << " [m]" << "\n";

      std::cout << std::setw(14) << "#B = " << Nb_ << "\n";
      std::cout << std::setw(14) << "D_b = " << Db_ << " [m^2*s^-1]" << "\n";
      std::cout << std::setw(14) << "r_b = " << rb_ << " [m]" << "\n";

      std::cout << std::setw(14) << "#C = " << Nc_ << "\n";
      std::cout << std::setw(14) << "D_c = " << Dc_ << " [m^2*s^-1]" << "\n";
      std::cout << std::setw(14) << "r_c = " << rc_ << " [m]" << "\n";

      std::cout << std::setw(14) << "ka = " << ka_ << " [m^3*s^-1]" << "\n";
      std::cout << std::setw(14) << "kd = " << kd_ << " [s^-1]" << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      sid_a = model_.add_species_type(SpeciesType("A", model_.get_def_structure_type_id(), Da_, ra_));
      sid_b = model_.add_species_type(SpeciesType("B", model_.get_def_structure_type_id(), Db_, rb_));
      sid_c = model_.add_species_type(SpeciesType("C", model_.get_def_structure_type_id(), Dc_, rc_));

      rules_.add_reaction_rule(ReactionRule(sid_a, sid_b, ka_, std::vector<SpeciesTypeID>{sid_c}));
      rules_.add_reaction_rule(ReactionRule(sid_c, kd_, std::vector<SpeciesTypeID>{sid_a, sid_b}));

      Simulation::SetupSimulation();

      world_.throwInParticles(sid_a, Na_, rng_);
      world_.throwInParticles(sid_b, Nb_, rng_);
      world_.throwInParticles(sid_c, Nc_, rng_);

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   int Na_ = 10;                                 // The number of particles of species type A
   double ra_ = 1e-9;                            // The radius of the A particles [m]
   double Da_ = 1e-12;                           // The diffusion constant D of the A particles [m^2*s^-1]

   int Nb_ = 10;                                 // The number of particles of species type B
   double rb_ = 1e-9;                            // The radius of the B particles [m]
   double Db_ = 1e-12;                           // The diffusion constant D of the B particles [m^2*s^-1]

   int Nc_ = 0;                                  // The number of particles of species type C.
   double rc_ = 1e-9;                            // The radius of the C particles [m]
   double Dc_ = 1e-12;                           // The diffusion constant D of the C particles [m^2*s^-1]

   double ka_ = 1e-16;                           // The association constant ka of the reaction [m^3*s^-1]
   double kd_ = 2e2;                             // The dissociation constant kd of the reaction [s^-1]

   // --------------------------------------------------------------------------------------------------------------------------------

   SpeciesTypeID sid_a;                          // The species type of the A particles
   SpeciesTypeID sid_b;                          // The species type of the B particles
   SpeciesTypeID sid_c;                          // The species type of the C particles

   // --------------------------------------------------------------------------------------------------------------------------------

};

#endif /* SIMULATION3P_HPP */
