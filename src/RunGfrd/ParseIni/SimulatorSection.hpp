#ifndef SIMULATOR_SECTION_HPP
#define SIMULATOR_SECTION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include "SectionBase.hpp"
#include "ParserExceptions.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

struct SimulatorSection final : SectionBase
{
   explicit SimulatorSection() : seed_(0), prepare_time_(0.0), end_time_(0.0), maintenance_step_(0), show_progress_(false) { }

   ~SimulatorSection() = default;

   // --------------------------------------------------------------------------------------------------------------------------------

   static std::string section_name() { return "Simulator"; }

   const std::string key_seed = "Seed";
   const std::string key_prepare_time = "PrepareTime";
   const std::string key_end_time = "EndTime";
   const std::string key_maintenance_step = "MaintenanceStep";
   const std::string key_maintenance_file = "MaintenanceFile";
   const std::string key_show_progress = "ShowProgress";

   int seed() const { return seed_; }
   double prepare_time() const { return prepare_time_; }
   double end_time() const { return end_time_; }
   size_t maintenance_step() const { return maintenance_step_; }
   std::string maintenance_file() const { return simstate_file_; }
   bool show_progress() const { return show_progress_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool set_keypair(const std::string& key, const std::string& value) override
   {
      if (key == key_seed) { seed_ = std::stoul(value, nullptr, 0); return true; }
      if (key == key_prepare_time) { prepare_time_ = std::stod(value); return true; }
      if (key == key_end_time) { end_time_ = std::stod(value); return true; }
      if (key == key_maintenance_step) { maintenance_step_ = std::stoi(value); return true; }
      if (key == key_maintenance_file) { simstate_file_ = value; return true; }
      if (key == key_show_progress) { show_progress_ = get_bool(value); return true; }
      THROW_EXCEPTION(illegal_section_key, "Key '" << key << "' not recognized.");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   int seed_;
   double prepare_time_;
   double end_time_;
   size_t maintenance_step_;
   std::string simstate_file_;
   bool show_progress_;
};

// --------------------------------------------------------------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& stream, const SimulatorSection& ss)
{
   stream << "[" << ss.section_name() << "]" << std::endl;
   if (ss.seed()) stream << ss.key_seed << " = 0x" << std::uppercase << std::hex << ss.seed() << std::dec << std::endl;
   if (ss.prepare_time()) stream << ss.key_prepare_time << " = " << ss.prepare_time() << std::endl;
   if (ss.end_time()) stream << ss.key_end_time << " = " << ss.end_time() << std::endl;
   if (ss.maintenance_step()) stream << ss.key_maintenance_step << " = " << ss.maintenance_step() << std::endl;
   if (ss.maintenance_step() && !ss.maintenance_file().empty()) stream << ss.key_maintenance_file << " = " << ss.maintenance_file() << std::endl;
   if (ss.show_progress()) stream << ss.key_show_progress << " = " << ss.show_progress() << std::endl;
   stream << std::endl;
   return stream;
}

#endif