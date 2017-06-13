#ifndef COPYNUMBERS_HPP
#define COPYNUMBERS_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <sstream>
#include <Logger.hpp>
#include "World.hpp"
#include "Event.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class CopyNumbers : public CustomAction
{
public:
   explicit CopyNumbers(const World& world, double interval = 1e-6) : CustomAction(interval), world_(world), log_(Log("CopyNummber"))
   {
      log_.set_stream(std::cout);
      log_.set_flags(Logger::logflags::None);
   }

   void do_action(double time) override
   {
      if (time == 0.0) print_header();
      print_row(time);
   }

   const char* type_name() const override { return "CopyNumbers"; }

   void set_output(std::ostream& stream) const { log_.set_stream(stream); }

private:
   void print_header() const
   {
      std::stringstream line;
      line << std::setw(12) << "Time";
      for (auto& s : world_.get_species())
         line << std::setw(8) << s.second.name();
      log_.info() << line.str();
   }

   void print_row(double time) const
   {
      std::stringstream line;
      line << std::scientific << std::setprecision(6) << std::setw(12) << time << std::fixed;
      for (auto& s : world_.get_species())
         line << std::setw(8) << world_.get_particle_ids(s.first).size();
      log_.info() << line.str();
   }

protected:
   const World& world_;
   Logger& log_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /* COPYNUMBERS_HPP */