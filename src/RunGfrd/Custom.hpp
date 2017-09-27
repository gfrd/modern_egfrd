#ifndef CUSTOM_HPP
#define CUSTOM_HPP

// --------------------------------------------------------------------------------------------------------------------------------

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
      return Simulation::PrepareSimulation();
   }

   void PrintSettings() override
   {
      Simulation::PrintSettings();
      std::cout << std::setw(14) << "input file = " << settingsfile_ << "\n";
      for (size_t i = 0; i < settings_.parameter_size(); i++)
         if (!settings_.get_parameter(i).empty()) std::cout << std::setw(14) << ("p" + std::to_string(i) + " = ") << settings_.get_parameter(i) << "\n";

      for (auto& species : settings_.getSpeciesTypeSections())
      {
         std::cout << std::setw(14) << "species = " << "'" << species.name() << "'" << ", D = " << species.D() << " [m^2*s^-1], r = " << species.r() << " [m]";
         if (species.v() != 0) std::cout << ", v = " << species.v();
         std::cout << "\n";
      }

      for (auto& rule : settings_.getReactionRuleSections())
      {
         std::cout << std::setw(14) << "rule = " << "'" << rule.rule() << "'" << ", k = " << rule.k() << "\n";
      }

      for (auto& pair : settings_.getParticlesSection().particles())
      {
         std::cout << std::setw(14) << "particle = " << "'" << pair.first << "'" << ", N = " << pair.second << "\n";
      }
   }

   bool SetupSimulation() override
   {
      set_species_section(settings_.getSpeciesTypeSections());
      set_reactionrule_section(settings_.getReactionRuleSections());
      
      Simulation::SetupSimulation();
      
      set_throw_particles_section(settings_.getParticlesSection());

      auto cns = settings_.getCopyNumbersSection();
      if (cns.mode() == CopyNumbersSection::modes::On)
         set_copynumbers_section(cns);
      auto pps = settings_.getParticlePositionsSection();
      if (pps.mode() == ParticlePositionSection::modes::On)
         set_particlepostions_section(pps);

      return true;
   }

   void PostPreSimulation() override
   {
      Simulation::PostPreSimulation();

      auto cns = settings_.getCopyNumbersSection();
      if (cns.mode() != CopyNumbersSection::modes::Off)
         set_copynumbers_section(cns);
      auto pps = settings_.getParticlePositionsSection();
      if (pps.mode() != ParticlePositionSection::modes::Off)
         set_particlepostions_section(pps);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   void set_copynumbers_section(const CopyNumbersSection& cns)
   {
      if (cns.type() == CopyNumbersSection::types::Average)
      {
         cna_ = std::make_unique<CopyNumbersAvg>(world_, rules_, cns.interval());
         simulator_->add_extrnal_event(0, cna_.get());
         simulator_->add_reaction_recorder(cna_.get());
         if (!cns.file().empty())
         {
            cfile_.open(cns.file(), std::fstream::in | std::fstream::out | std::fstream::trunc);
            cna_->set_output(cfile_);
         }
      }
      else
      {
         cni_ = std::make_unique<CopyNumbersInst>(world_, cns.interval());
         simulator_->add_extrnal_event(0, cni_.get());
         if (!cns.file().empty())
         {
            cfile_.open(cns.file(), std::fstream::in | std::fstream::out | std::fstream::trunc);
            cni_->set_output(cfile_);
         }
      }
   }

   void set_particlepostions_section(const ParticlePositionSection& pps)
   {
      pp_ = std::make_unique<ParticlePositions>(simulator_, pps.interval());
      simulator_->add_extrnal_event(0, pp_.get());
      if (!pps.file().empty())
      {
         pfile_.open(pps.file(), std::fstream::in | std::fstream::out | std::fstream::trunc);
         pp_->set_output(pfile_);
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
   std::fstream pfile_;
   std::unique_ptr<CopyNumbersInst> cni_;
   std::unique_ptr<CopyNumbersAvg> cna_;
   std::unique_ptr<ParticlePositions> pp_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /*CUSTOM_HPP*/