#ifndef POWER_SPECTRUM_HPP
#define POWER_SPECTRUM_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <fstream>

// --------------------------------------------------------------------------------------------------------------------------------

using state_time_pair_arr = std::vector<std::pair<bool, double>>;
using vector2_arr = std::vector<Vector2>;
using double_arr = std::vector<double>;

// --------------------------------------------------------------------------------------------------------------------------------

class ReactionRecorder : public reaction_recorder_log
{
public:
   explicit ReactionRecorder(const ReactionRuleCollection& reaction_rules, const Model& model) : reaction_recorder_log(reaction_rules, model), count_(0), time_(0) { }

   size_t count() const { return count_; }
   double final_time() const { return time_; }

protected:
   void StoreReaction(double time, ReactionRuleID rid, ParticleID r1, ParticleID r2, ParticleID p1, ParticleID p2) override
   {
      reaction_recorder_log::StoreReaction(time, rid, r1, r2, p1, p2);
      count_++;
      time_ = time;
   }

private:
   size_t count_;
   double time_;
};


// --------------------------------------------------------------------------------------------------------------------------------

// Calculate power spectrum of bound/unbound particle.
class PowerSpectrum : public Simulation3P
{
   const int averages_count_ = 1 << 6;
   const int distribution_count_ = 1 << 15;

public:

   explicit PowerSpectrum() : Simulation3P(), rfilename_("reactionrecord.txt")
   {
      end_time_ = 60;
      Na_ = 1; Nb_ = 241; Nc_ = 0;
      ka_ = 9.16639e-19; kd_ = 220.8;
      world_size_ = 1e-6;
      ra_ = rb_ = 5e-9; rc_ = 1e-8;
      Da_ = 0; Db_ = 1e-12; Dc_ = 0;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "PowerSpectrum"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "rr" && args.isvalue(i + 1)) rfilename_ = args.option(++i);
      else return Simulation3P::HandleCommandArguments(i, args);
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      Simulation3P::print_usage();
      std::cout << "      -rr filename     Output file for ReactionRecord (default 'reactionrecord.txt')\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   void PrintSettings() override
   {
      Simulation3P::PrintSettings();
      std::cout << std::setw(15) << "react.log = '" << rfilename_ << "'\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      Simulation3P::SetupSimulation();
      rrec_ = std::make_unique<ReactionRecorder>(rules_, model_);
      simulator_->add_reaction_recorder(rrec_.get());
      rfile_.open(rfilename_, std::fstream::in | std::fstream::out | std::fstream::trunc);
      rrec_->set_output(rfile_);

      // check settings!
      if (Na_ != 1 || Nc_ != 0) { Log("RunGfrd").fatal() << "PowerSpectrum can only be calculated when Na = 1 and Nc = 0."; failed_ = true; }

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PostSimulation() override
   {
      if (rrec_->count() > 0) SavePowerSpectrum();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   void SavePowerSpectrum()
   {
      auto omega_max = 10.0;
      auto final_time = rrec_->final_time();
      auto omega_min = std::log10(2.0 * M_PI * omega_max / final_time);

      Log("RunGfrd").info() << "Calculating power spectrum";
      auto omega = logDistribution(omega_min, omega_max);
      auto fourier = fourierTransform(omega);
      auto power_spectrum = powerSpectrum(fourier, final_time);

      std::cout << "\n\nomega [1/s]\tP[s]\n";
      for (int i = 0; i < distribution_count_; ++i)
         std::cout << omega[i] << "\t" << power_spectrum[i] << "\n";

      Log("RunGfrd").info() << "Calculating averages";
      auto omega_average = averageDistribution(omega_min, omega_max);
      auto unweigthed_average = unweightedAverage(power_spectrum);
      auto weighted_average = weightedAverage(power_spectrum, omega);

      std::cout << "\n\nomega [1/s]\tP_unweighted[s]\tP_weigthed[s]\n";
      for (int i = 0; i < distribution_count_ / averages_count_; ++i)
         std::cout << omega_average[i] << "\t" << unweigthed_average[i] << "\t" << weighted_average[i] << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double_arr logDistribution(double omega_min, double omega_max) const
   {
      double_arr omega(distribution_count_);
      auto weight = (omega_max - omega_min) / distribution_count_;
      for (auto index = 0; index < distribution_count_; ++index)
         omega[index] = std::pow(10, (index + 1.0) * weight + omega_min);

      return omega;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double_arr averageDistribution(double omega_min, double omega_max) const
   {
      auto size = static_cast<int>(distribution_count_ / averages_count_);
      double_arr omega_average(size);

      auto weight = (omega_max - omega_min) / distribution_count_ * averages_count_;
      for (int index = 0; index < size; ++index)
         omega_average[index] = std::pow(10, (index + 1.0) * weight + omega_min);

      return omega_average;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   vector2_arr fourierTransform(const double_arr& omega)
   {
      auto states_size = rrec_->count();
      auto omega_size = omega.size();
      vector2_arr fourier(omega_size);

      rfile_.seekp(0); // parse the reaction record
      std::string line;
      while (line.find("Time [s]") == std::string::npos)
         std::getline(rfile_, line);

      auto print = std::max(10ull, states_size / 100ull);
      for (size_t i = 0; i < states_size; ++i)
      {
         if (i && i%print == 0) std::cerr << "\r" << (i * 100ull / states_size) << "% ";

         std::getline(rfile_, line);
         std::istringstream in(line);
         double time;
         idtype rid;
         in >> std::scientific >> time >> rid;

         double flip_state = rid == 1 ? 1 : -1;    // bind +1, unbind -1
         for (size_t j = 0; j < omega_size; ++j)
         {
            auto arg = omega[j] * time;
            auto coef = flip_state / omega[j];
            fourier[j] += coef * Vector2(-std::sin(arg), std::cos(arg));
         }
      }

      std::cerr << "\n";
      return fourier;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   static double_arr powerSpectrum(const vector2_arr& fourier, double max_time)
   {
      double_arr psd;
      psd.reserve(fourier.size());
      for (const auto& i : fourier)
         psd.emplace_back(Vector2::dot(i, i) / max_time);

      return psd;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double_arr unweightedAverage(const double_arr& spectrum) const
   {
      size_t size = distribution_count_ / averages_count_;
      double_arr average(size);

      average.resize(size);
      for (size_t index = 0; index < size; index++)
      {
         auto sum = 0.0;
         for (auto lag = 0; lag < averages_count_; lag++)
            sum += spectrum[averages_count_ * index + lag];

         average[index] = sum / averages_count_;
      }

      return average;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double_arr weightedAverage(const double_arr& spectrum, const double_arr& omega) const
   {
      auto size = distribution_count_ / averages_count_;
      double_arr average(size);

      for (auto index = 0; index < size; index++)
      {
         auto sum = 0.0;
         for (auto lag = 0; lag < averages_count_ - 1; lag++)
         {
            auto numerator = spectrum[averages_count_*index + lag] + spectrum[averages_count_ * index + lag + 1];
            auto denominator = 2.0 * (omega[averages_count_ * index + lag + 1] - omega[averages_count_ * index + lag]);
            sum += numerator / denominator;
         }
         average[index] = sum / (omega[(index + 1) * averages_count_ - 1] - omega[index * averages_count_]);
      }

      return average;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string rfilename_;
   std::fstream rfile_;
   std::unique_ptr<ReactionRecorder> rrec_;

   // --------------------------------------------------------------------------------------------------------------------------------
};

#endif /* POWER_SPECTRUM_HPP */
