#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <Logger.hpp>
#include "Event.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Progress : public CustomAction
{
public:
   Progress(const std::unique_ptr<EGFRDSimulator>& simulator, double end_time, int length = 40) :
   CustomAction(end_time / (end_time > 0.1 ? 1000.0 : 200.0)), simulator_(simulator), log_(Log("Progress")), end_time_(end_time), length_(length)
   {
      log_.set_stream(std::clog);
      log_.set_flags(Logger::logflags::Name | Logger::logflags::Separator);
      log_.set_postfix("\r");
   }

   void do_action(double time) override
   {
      print_progress(time);
   }

   const char* type_name() const override { return "Progress"; }

private:
   void print_progress(double time) const
   {
        const double progress = time / end_time_;
        const auto stars = static_cast<int>(std::round(length_ * progress));

        std::array<uint, 3> domain_type_count{ 0,0,0 };
        for (auto& domain : simulator_->get_domains())
        {
           domain_type_count[static_cast<int>(domain.second->multiplicity()) - 1]++;
        }
      
        std::stringstream line;
        line << "[";
        for (int i = 0; i < stars; i++) line << "*";
        for (int i = stars; i < length_; i++) line << ".";
        line << "] " << std::fixed << std::setprecision(1) << progress * 100.0 << " %\t";
        line << "(steps: " << simulator_->num_steps() << std::setprecision(12) << ", dt per step: " << double(simulator_->time()) / simulator_->num_steps() << ")\t";
        line << "{S: " << domain_type_count[0] << " P: " << domain_type_count[1] << " M: " << domain_type_count[2] << "}";
        log_.info() << line.str();
   }

   Logger log_;
   double end_time_;
   int length_;
   const std::unique_ptr<EGFRDSimulator>& simulator_;
};

// --------------------------------------------------------------------------------------------------------------------------------
