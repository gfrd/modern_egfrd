#ifndef CUSTOM_HPP
#define CUSTOM_HPP
#include "ParseIni/SimulatorSettings.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class Custom : public Simulation
{
public:

   explicit Custom() noexcept : settingsfile_() {}

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Custom"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "f" && args.isvalue(i + 1)) settingsfile_ = args.option(++i);
      else if (args.isparam(i) && *args.option(i).begin() == 'p' && args.isvalue(i + 1))
      {
         size_t p = args.option(i).length() == 2 ? args.option(i)[1] - '0' : -1;
         if (p < 0 || p >= settings_.parameter_size()) return true;
         settings_.set_parameter(p, args.option(++i));
      }
      else return true;
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      std::cout << "      -f filename      Settings/Model INI file\n";
      std::cout << "      -p<0..9> value   Give parameter N a value (use <paramN> in INI file)\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   bool PrepareSimulation() override
   {
      if (settingsfile_.empty()) { Log("RunGfrd").fatal() << "No input file! Specify with -f command option"; return false; }

      std::ifstream stream(settingsfile_);
      if (!stream.is_open()) { Log("RunGfrd").fatal() << "Failed to load input file!"; return false; }
      stream >> settings_;
      stream.close();

      set_simulator_section(settings_.getSimulatorSection());
      set_world_section(settings_.getWorldSection());
      set_species_section(settings_.getSpeciesTypeSections());
      set_reactionrule_section(settings_.getReactionRuleSections());
      return Simulation::PrepareSimulation();
   }

   void PrintSettings() override
   {
      Simulation::PrintSettings();
      std::cout << std::setw(14) << "input file = " << settingsfile_ << "\n";
      for (size_t i=0; i< settings_.parameter_size(); i++)
         if (!settings_.get_parameter(i).empty()) std::cout << std::setw(14) << ("p" + std::to_string(i) + " = ") << settings_.get_parameter(i) << "\n";
   }

   bool SetupSimulation() override
   {
      Simulation::SetupSimulation();
      set_throw_particles_section(settings_.getParticlesSection());

      auto cns = settings_.getCopyNumbersSection();
      if (cns.mode() == CopyNumbersSection::modes::On)
         set_copynumbers_section(cns);

      return true;
   }

   void PostPreSimulation() override
   {
      Simulation::PostPreSimulation();

      auto cns = settings_.getCopyNumbersSection();
      if (cns.mode() != CopyNumbersSection::modes::Off)
         set_copynumbers_section(cns);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   void set_copynumbers_section(const CopyNumbersSection& cns)
   {
      cn_ = std::make_unique<CopyNumbers>(world_, cns.interval());
      simulator_->add_extrnal_event(0, cn_.get());
      if (!cns.file().empty())
      {
         cfile_.open(cns.file(), std::fstream::in | std::fstream::out | std::fstream::trunc);
         cn_->set_output(cfile_);
      }
   }

   void set_simulator_section(const SimulatorSection& section)
   {
      if (section.seed()) rng_.seed(section.seed());
      prep_time_ = section.prepare_time();
      end_time_ = section.end_time();
      show_progress_ = section.show_progress();
      maintenance_step_ = section.maintenance_step();
      simstate_file_ = section.maintenance_file();
   }

   void set_species_section(const std::vector<SpeciesTypeSection>& sections)
   {
      for (auto& section : sections)
      {
         auto species = section.create_species(model_.get_def_structure_type_id());
         model_.add_species_type(species);
      }
   }

   void set_world_section(const WorldSection& section)
   {
      world_size_ = section.world_size();
      THROW_UNLESS_MSG(illegal_matrix_size, world_.matrix_size()[0] == section.matrix_space(), "invalid matrix space x-size!");
      THROW_UNLESS_MSG(illegal_matrix_size, world_.matrix_size()[1] == section.matrix_space(), "invalid matrix space y-size!");
      THROW_UNLESS_MSG(illegal_matrix_size, world_.matrix_size()[2] == section.matrix_space(), "invalid matrix space z-size!");
   }

   void set_reactionrule_section(const std::vector<ReactionRuleSection>& sections)
   {
      for (auto& section : sections)
         rules_.add_reaction_rule(section.create_reaction_rule(model_));
   }

   void set_throw_particles_section(const ParticlesSection& section)
   {
      section.add_particles_to_world(model_, world_, rng_);
   }

   // --------------------------------------------------------------------------------------------------------------------------------
   
   std::string settingsfile_;
   SimulatorSettings settings_;

   std::fstream cfile_;
   std::unique_ptr<CopyNumbers> cn_;

};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /*CUSTOM_HPP*/