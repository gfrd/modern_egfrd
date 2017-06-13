#ifndef EQUILIBRIUM_HPP
#define EQUILIBRIUM_HPP

#include <CopyNumbers.hpp>
#include "Simulation3P.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class CopyNumbersAvg : public CopyNumbers
{
public:
   CopyNumbersAvg(const World& w, double interval, SpeciesTypeID sid) : CopyNumbers(w, interval), sid_(sid), count_(0), sum_(0) { }

   void do_action(double time) override
   {
      CopyNumbers::do_action(time);
      sum_ += world_.get_particle_ids(sid_).size();
      count_++;
   }

   void print_avg() const
   {
      if (count_ > 0) log_.info() << "Average NC  = " << static_cast<double>(sum_) / count_;
      else log_.info() << "Average NC  = N.A. ";
   }

   const char* type_name() const override { return "CopyNumbersAvg"; }

private:
   SpeciesTypeID sid_;     // Calculate average number of particles with id
   size_t count_, sum_;
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
      cn_ = std::make_unique<CopyNumbersAvg>(world_, ctime_, sid_c);
      simulator_->add_extrnal_event(0, cn_.get());
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
      if (cn_.get()) cn_->print_avg();
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
   std::unique_ptr<CopyNumbersAvg> cn_;
   double ctime_ = 1E-1;


   // --------------------------------------------------------------------------------------------------------------------------------

};

#endif /* EQUILIBRIUM_HPP */