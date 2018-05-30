#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <Logger.hpp>
#include "Event.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Progress : public CustomAction
{
public:
   Progress(double end_time, int length = 40) : CustomAction(end_time / (end_time > 0.1 ? 1000.0 : 200.0)), log_(Log("Progress")), end_time_(end_time), length_(length)
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
      std::stringstream line;
      line << "[";
      for (int i = 0; i < stars; i++) line << "*";
      for (int i = stars; i < length_; i++) line << ".";
      line << "] " << std::fixed << std::setprecision(1) << progress * 100.0 << " %";
      log_.info() << line.str();
   }

   Logger log_;
   double end_time_;
   int length_;
};

// --------------------------------------------------------------------------------------------------------------------------------
