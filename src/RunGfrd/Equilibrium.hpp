#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include <CopyNumbers.hpp>
#include "Simulation3P.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class CopyNumberFull : public CopyNumbersAvg
{
   // Computes the average number of particles of type 'sid' over the whole simulation time

public:
   CopyNumberFull(const World& w, const ReactionRuleCollection& rules, double interval, SpeciesTypeID sid) : CopyNumbersAvg(w, rules, interval), sid_(sid) { }

   void do_action(double time) override
   {
      CopyNumbersAvg::do_action(time);
      if (time == 0) accu_sid_.reset(world_.get_particle_ids(sid_).size());
   }

   void print_avg(double time)
   {
      auto avg = accu_sid_.update_time(time);
      if (time > 0) log_.info() << "Average NC  = " << avg;
      else log_.info() << "Averag2 NC  = N.A. ";
   }

   const char* type_name() const override { return "CopyNumberFull"; }

protected:
   void StoreReaction(double time, ReactionRuleID rid, ParticleID r1, ParticleID r2, ParticleID p1, ParticleID p2) override
   {
      CopyNumbersAvg::StoreReaction(time, rid, r1, r2, p1, p2);

      size_t count = world_.get_particle_ids(sid_).size();
      accu_sid_.update_count(time, count);
   }

private:
   SpeciesTypeID sid_;     // Calculate average number of particles with id
   avg_accu accu_sid_;
};

// ================================================================================================================================

class Equilibrium : public Simulation3P
{
public:
   Equilibrium() : Simulation3P()
   {
      Na_ = 100;
      Nb_ = 100;
      Nc_ = 0;
      world_size_ = 3.42e-6;        // 40 femto Liter = 4E-17 m3
      equil_start_ = false;
      prep_time_ = 0.0;
      ka_ = 1e-19;                  // The association constant ka of the reaction [m^3*s^-1]
      kd_ = 2e-2;                   // The dissociation constant kd of the reaction [s^-1]
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Equilibrium"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "cn" && args.isvalue(i + 1)) cfilename_ = args.option(++i);
      else if (args.isparam(i) && args.option(i) == "ci" && args.isvalue(i + 1)) ctime_ = args.value(++i);
      else return Simulation3P::HandleCommandArguments(i, args);
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      Simulation3P::print_usage();
      std::cout << "      -cn filename      Output file for CopyNumbers (default 'std::out')\n";
      std::cout << "      -ci x.xx          CopyNumbers output interval (default 1E-3 sec)\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------
protected:

   bool PrepareSimulation() override
   {
      equil_start_ = Nc_ < 0;
      if (equil_start_)      // start in equilibrium (with given Na and Nb, where Nc = 0, calculate new Na,Nb,Nc, total particles stays the same plus/minus one or two )
      {
         Nc_ = 0;
         int Na = Na_;
         while (true)
         {
            Nc_ = static_cast<int>(std::round(AverageParticles()));
            double delta = Na - (Na_ + Nc_);
            if (std::fabs(delta) < 0.66667) break;
            int d = static_cast<int>(std::round(0.75*delta));  Na_ += d;  Nb_ += d;
         }
      }
      return Simulation::PrepareSimulation();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() override
   {
      Simulation3P::PrintSettings();
      std::cout << std::setw(14) << "equalized = " << (equil_start_ ? "yes" : "no") << "\n";
      std::cout << std::setw(14) << "<NC> = " << AverageParticles() << "\n";
      std::cout << std::setw(15) << "copy.num. = '" << (cfilename_.empty() ? "std::out" : cfilename_) << "'\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void AddCopyNumbers()
   {
      cn_ = std::make_unique<CopyNumberFull>(world_, rules_, ctime_, sid_c);
      simulator_->add_extrnal_event(0, cn_.get());
      simulator_->add_reaction_recorder(cn_.get());

      if (!cfilename_.empty())
      {
         cfile_.open(cfilename_, std::fstream::in | std::fstream::out | std::fstream::trunc);
         cn_->set_output(cfile_);
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      Simulation3P::SetupSimulation();

      if (prep_time_ == 0)    // without pre-sim phase, we need to init the copynumber here.
         AddCopyNumbers();
      return true;
   }
   
   // --------------------------------------------------------------------------------------------------------------------------------

   void PostPreSimulation() override
   {
      Simulation3P::PostPreSimulation();
      AddCopyNumbers();           // with pre-sim phase, init the copynumber for main-sim phase.
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PostSimulation() override
   {
      if (cn_.get()) cn_->print_avg(simulator_->time());
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   
   // Calculates the analytic average number of C particles in equilibrium.
   double AverageParticles() const
   {
      int Na = Na_ + Nc_;     // disassociated c particles into a and b
      int Nb = Nb_ + Nc_;

      double kdv = kd_ / ka_ * std::pow(world_size_, 3);
      double sum = Na + Nb + kdv;
      return 0.5 * (sum - std::sqrt(sum * sum - 4 * Na * Nb));
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool equil_start_;
   std::string cfilename_;
   std::fstream cfile_;
   std::unique_ptr<CopyNumberFull> cn_;
   double ctime_ = 1E-1;


   // --------------------------------------------------------------------------------------------------------------------------------

};

#endif /* EQUILIBRIUM_HPP */