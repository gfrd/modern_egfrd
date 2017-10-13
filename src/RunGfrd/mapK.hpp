#ifndef MAPK_HPP
#define MAPK_HPP

#include <CopyNumbers.hpp>
#include "Simulation.hpp"
#include "convert.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class MapK : public Simulation
{
public:
   explicit MapK() : Simulation(), cfilename_(), rfilename_("reactionrecord.txt")
   {
      prep_time_ = 0.1;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "MapK"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "D" && args.isvalue(i + 1)) D_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "Nr" && args.isvalue(i + 1)) Nratio_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "Trel" && args.isvalue(i + 1)) Trel_ = args.value(++i);
      else if (args.isparam(i) && args.option(i) == "rr" && args.isvalue(i + 1)) rfilename_ = args.option(++i);
      else if (args.isparam(i) && args.option(i) == "cn" && args.isvalue(i + 1)) cfilename_ = args.option(++i);
      else if (args.isparam(i) && args.option(i) == "ci" && args.isvalue(i + 1)) ctime_ = args.value(++i);
      else return Simulation::HandleCommandArguments(i, args);
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      Simulation::print_usage();
      std::cout << "      -D x.xx          Diffusion coefficient of all particles in  m^2*s^-1 (default 1E-12)\n";
      std::cout << "      -Nr x.xx         Ratio of dephosporylated and double phosphorylated proteins (default 0.5)\n";
      std::cout << "      -Trel x.xx       Relaxation time of the KK* and P* enzymes in seconds (default 1E-2)\n";
      std::cout << "      -rr filename      Output file for ReactionRecord (default 'reactionrecord.txt')\n";
      std::cout << "      -cn filename      Output file for CopyNumbers (default 'std::out')\n";
      std::cout << "      -ci x.xx          CopyNumbers output interval (default 0.1 sec)\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   bool PrepareSimulation() override
   {
      world_size_ = std::pow(volume_ * 1e-3, 1.0 / 3.0);

      N_ = static_cast<uint>(convert::particles_in_volume_liter(conS_, volume_));
      N_Kpp = static_cast<uint>(N_ * Nratio_);
      N_K = N_ - N_Kpp;
      N_KK = static_cast<uint>(convert::particles_in_volume_liter(0.5 * conE_, volume_));
      N_P = static_cast<uint>(convert::particles_in_volume_liter(0.5 * conE_, volume_));

      k7 = std::log(2.0) / Trel_;

      return Simulation::PrepareSimulation();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() override
   {
      Simulation::PrintSettings();
      std::cout << std::setw(14) << "volume = " << std::setprecision(6) << volume_ << " [L]\n";
      std::cout << std::setw(14) << "radius = " << radius_ << " [m]\n";
      std::cout << std::setw(14) << "con Ez = " << conE_ << " [M]\n";
      std::cout << std::setw(14) << "con SS = " << conS_ << " [M]\n";
      std::cout << std::setw(14) << "k1 = " << k1_ << " [nM^-1*s^-1]\n";
      std::cout << std::setw(14) << "k2 = " << k2_ << " [s^-1]\n";
      std::cout << std::setw(14) << "k3 = " << k3_ << " [s^-1]\n";
      std::cout << std::setw(14) << "k4 = " << k4_ << " [nM^-1*s^-1]\n";
      std::cout << std::setw(14) << "k5 = " << k5_ << " [s^-1]\n";
      std::cout << std::setw(14) << "k6 = " << k6_ << " [s^-1]\n";

      std::cout << std::setw(14) << "D = " << D_ << " [m^2*s^-1]\n";
      std::cout << std::setw(14) << "N ratio = " << Nratio_ << "\n";
      std::cout << std::setw(14) << "Trel = " << std::setprecision(6) << Trel_ << " [s]\n";
      std::cout << std::setw(14) << "k7 = " << k7 << " [s^-1]\n";

      std::cout << std::setw(14) << "N_K   = " << N_K << std::endl;
      std::cout << std::setw(14) << "N_Kpp = " << N_Kpp << std::endl;
      std::cout << std::setw(14) << "N_KK  = " << N_KK << std::endl;
      std::cout << std::setw(14) << "N_P   = " << N_P << std::endl;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      StructureTypeID wsid = model_.get_def_structure_type_id();

      sK = model_.add_species_type(SpeciesType("K", wsid, D_, radius_));
      sKK_K = model_.add_species_type(SpeciesType("K_KK", wsid, D_, radius_));
      sKp = model_.add_species_type(SpeciesType("Kp", wsid, D_, radius_));
      sKK_Kp = model_.add_species_type(SpeciesType("Kp_KK", wsid, D_, radius_));
      sP_Kp = model_.add_species_type(SpeciesType("Kp_P", wsid, D_, radius_));
      sKpp = model_.add_species_type(SpeciesType("Kpp", wsid, D_, radius_));
      sP_Kpp = model_.add_species_type(SpeciesType("Kpp_P", wsid, D_, radius_));

      sKK = model_.add_species_type(SpeciesType("KK", wsid, D_, radius_));
      sKKi = model_.add_species_type(SpeciesType("KKi", wsid, D_, radius_));
      sP = model_.add_species_type(SpeciesType("P", wsid, D_, radius_));
      sPi = model_.add_species_type(SpeciesType("Pi", wsid, D_, radius_));

      Simulation::SetupSimulation();

      world_.throwInParticles(sK, N_K, rng_);
      world_.throwInParticles(sKpp, N_Kpp, rng_);
      world_.throwInParticles(sKK, N_KK, rng_);
      world_.throwInParticles(sP, N_P, rng_);

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PostPreSimulation() override
   {
      Simulation::PostPreSimulation();
      auto kk1 = convert::per_nM_per_sec_to_m3_per_sec(k1_);
      auto kk4 = convert::per_nM_per_sec_to_m3_per_sec(k4_);

      //[1]
      rules_.add_reaction_rule(ReactionRule(sKK, sK, kk1, std::vector<SpeciesTypeID>{ sKK_K }));
      rules_.add_reaction_rule(ReactionRule(sKK_K, k2_, std::vector<SpeciesTypeID>{ sKK, sK }));
      rules_.add_reaction_rule(ReactionRule(sKK_K, k3_, std::vector<SpeciesTypeID>{ sKKi, sKp }));

      //[2]
      rules_.add_reaction_rule(ReactionRule(sKK, sKp, kk4, std::vector<SpeciesTypeID>{ sKK_Kp }));
      rules_.add_reaction_rule(ReactionRule(sKK_Kp, k5_, std::vector<SpeciesTypeID>{ sKK, sKp }));
      rules_.add_reaction_rule(ReactionRule(sKK_Kp, k6_, std::vector<SpeciesTypeID>{ sKKi, sKpp }));

      //[3]
      rules_.add_reaction_rule(ReactionRule(sP, sKpp, kk1, std::vector<SpeciesTypeID>{ sP_Kpp }));
      rules_.add_reaction_rule(ReactionRule(sP_Kpp, k2_, std::vector<SpeciesTypeID>{ sP, sKpp }));
      rules_.add_reaction_rule(ReactionRule(sP_Kpp, k3_, std::vector<SpeciesTypeID>{ sPi, sKp }));

      //[4]
      rules_.add_reaction_rule(ReactionRule(sP, sKp, kk4, std::vector<SpeciesTypeID>{ sP_Kp }));
      rules_.add_reaction_rule(ReactionRule(sP_Kp, k5_, std::vector<SpeciesTypeID>{  sP, sKp }));
      rules_.add_reaction_rule(ReactionRule(sP_Kp, k6_, std::vector<SpeciesTypeID>{ sPi, sK }));

      //[5]
      rules_.add_reaction_rule(ReactionRule(sKKi, k7, std::vector<SpeciesTypeID>{ sKK }));
      rules_.add_reaction_rule(ReactionRule(sPi, k7, std::vector<SpeciesTypeID>{ sP }));

      // copy nunmbers
      cn_ = std::make_unique<CopyNumbersInst>(world_, ctime_);
      simulator_->add_extrnal_event(0, cn_.get());
      if (!cfilename_.empty())
      {
         cfile_.open(cfilename_, std::fstream::in | std::fstream::out | std::fstream::trunc);
         cn_->set_output(cfile_);
      }

      // reaction record
      rrec_ = std::make_unique<reaction_recorder_log>(rules_, model_);
      simulator_->add_reaction_recorder(rrec_.get());
      rfile_.open(rfilename_, std::fstream::in | std::fstream::out | std::fstream::trunc);
      rrec_->set_output(rfile_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   // Model const parameters
   const double volume_ = 1E-15;                     // one femtolitre or a cubic micro liter
   const double radius_ = 2.5e-9;                    // particle radius
   const double conE_ = 100e-9;                      // Enzyme concentration (total KK + P)
   const double conS_ = 200e-9;                      // Substrate concentration (total K + Kp + Kpp)
   const double k1_ = 0.027;                         // intrinsic rate in nM^-1.s^-1
   const double k2_ = 1.35;                          // intrinsic rate in s^-1
   const double k3_ = 1.5;                           // intrinsic rate in s^-1
   const double k4_ = 0.056;                         // intrinsic rate in nM^-1.s^-1
   const double k5_ = 1.73;                          // intrinsic rate in s^-1
   const double k6_ = 15.0;                          // intrinsic rate in s^-1

   // Model input parameters
   double Trel_ = 10E-3;                              // Enzyme reactivation time in seconds
   double D_ = 1E-12;                                 // Diffusion constant (all particles) in m^2/s
   double Nratio_ = 0.5;                              // Start ratio between K and Kpp

   // Model calculated parameters
   double k7;
   uint N_, N_Kpp, N_K, N_KK, N_P;

   double ctime_ = 0.1;
   std::string cfilename_;
   std::fstream cfile_;
   std::unique_ptr<CopyNumbersInst> cn_;
   std::string rfilename_;
   std::fstream rfile_;
   std::unique_ptr<reaction_recorder_log> rrec_;

   SpeciesTypeID sK, sKp, sKpp;
   SpeciesTypeID sKK, sKKi;
   SpeciesTypeID sP, sPi;
   SpeciesTypeID sKK_K, sKK_Kp, sP_Kpp, sP_Kp;

   // --------------------------------------------------------------------------------------------------------------------------------

};




// --------------------------------------------------------------------------------------------------------------------------------

#endif /* MAPK_HPP */