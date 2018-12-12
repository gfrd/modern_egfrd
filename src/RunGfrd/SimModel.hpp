#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

#include <utility>
#include "SimulatorSettings.hpp"
#include "CopyNumbers.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

class SimModel : public Simulation
{
public:

   // --------------------------------------------------------------------------------------------------------------------------------

   explicit SimModel(std::string file) noexcept : file_(std::move(file)) {}

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Model"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   int HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "d" && args.isvalue(i + 1)) settings_.add_variable(args.option(++i));
      else return Simulation::HandleCommandArguments(i, args);
      return -1;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      std::cout << "  RunGfrd <path-to-model> [ options ]" << std::endl;
      std::cout << "        [-h,-?,--help  ]      Print command line usage information" << std::endl;
      std::cout << "        [-d var=value  ]      Define/Overwrite model value variable\n";
      std::cout << "        [-d var=$'text']      Define/Overwrite model string variable\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   bool PrepareSimulation() override
   {
      std::ifstream stream(file_);
      if (!stream.is_open()) { Log("RunGfrd").fatal() << "Failed to load input file!"; return false; }
      stream >> settings_;
      stream.close();

      {
         auto &section = settings_.getSimulatorSection();
         if (section.seed()) rng_.seed(seed_ = section.seed());
         prep_time_ = section.prepare_time();
         end_time_ = section.end_time();
         maintenance_step_ = section.maintenance_step();
         simstate_file_ = section.maintenance_file();
      }

      {
         auto& section(settings_.getWorldSection());
         world_size_ = section.world_size();
         THROW_UNLESS_MSG(illegal_size, world_size_ > 0, "invalid world-size!");
         THROW_UNLESS_MSG(illegal_size, world_.matrix_size()[0] == section.matrix_space(), "invalid matrix space x-size!");
         THROW_UNLESS_MSG(illegal_size, world_.matrix_size()[1] == section.matrix_space(), "invalid matrix space y-size!");
         THROW_UNLESS_MSG(illegal_size, world_.matrix_size()[2] == section.matrix_space(), "invalid matrix space z-size!");
      }

      return Simulation::PrepareSimulation();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PrintSettings() override
   {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      auto time_local = std::localtime(&time);
      std::cout << std::setw(14) << "time local = " << std::asctime(time_local);
      std::cout << std::setw(14) << "model file = " << file_ << "\n";

      std::cout << "\n";

      settings_.getVariablesSection().PrintSettings();

      Simulation::PrintSettings();

      std::cout << "\n";

      for (auto& species : settings_.getSpeciesTypeSections())
         species.PrintSettings();

      std::cout << "\n";

      for (auto& rule : settings_.getReactionRuleSections())
         rule.PrintSettings();

      std::cout << "\n";

      for (auto& section : settings_.getParticlesSections())
         section.PrintSettings();

      std::cout << "\n";

      auto cns = settings_.getCopyNumbersSection();
      if (cns != nullptr) cns->PrintSettings();
      auto pps = settings_.getParticlePositionsSection();
      if (pps != nullptr) pps->PrintSettings();
      auto rrs = settings_.getReactionRecordSection();
      if (rrs != nullptr) rrs->PrintSettings();
      auto ps = settings_.getProgressSection();
      if (ps != nullptr) ps->PrintSettings();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      auto& vars = settings_.getVariablesSection();

      {
         auto& sections = settings_.getSpeciesTypeSections();
         for (auto& section : sections)
            section.create_species(model_, vars);
      }

      Simulation::SetupSimulation();
      set_reactionrule_section(true);
      add_particles_section(true);

      auto cns = settings_.getCopyNumbersSection();
      if (cns != nullptr && test_mode_setup(true, cns->mode())) set_copynumbers_section(*cns);
      auto pps = settings_.getParticlePositionsSection();
      if (pps != nullptr && test_mode_setup(true, pps->mode())) set_particlepostions_section(*pps);
      auto rrs = settings_.getReactionRecordSection();
      if (rrs != nullptr && test_mode_setup(true, rrs->mode())) set_reactionrecord_section(*rrs);
      auto ps = settings_.getProgressSection();
      if (ps != nullptr && test_mode_setup(true, ps->mode())) set_progress_section(true, ps->column_width());

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void PostPreSimulation() override
   {
      Simulation::PostPreSimulation();
      set_reactionrule_section(false);
      add_particles_section(false);

      auto cns = settings_.getCopyNumbersSection();
      if (cns != nullptr && test_mode_setup(false, cns->mode(), true)) set_copynumbers_section(*cns);
      auto pps = settings_.getParticlePositionsSection();
      if (pps != nullptr && test_mode_setup(false, pps->mode(), true)) set_particlepostions_section(*pps);
      auto rrs = settings_.getReactionRecordSection();
      if (rrs != nullptr && test_mode_setup(false, rrs->mode(), true)) set_reactionrecord_section(*rrs);
      auto ps = settings_.getProgressSection();
      if (ps != nullptr && test_mode_setup(false, ps->mode(), true)) set_progress_section(false, ps->column_width());
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:

   bool test_mode_setup(bool first, SectionModeBase::modes mode, bool both = false) const
   {
      switch (mode)
      {
      case SectionModeBase::modes::Off: return false;
      case SectionModeBase::modes::On: return first || both;
      case SectionModeBase::modes::Run: return prepare_time() == 0 || !first;
      }
      THROW_EXCEPTION(illegal_state, "unknown mode");
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void set_copynumbers_section(const CopyNumbersSection& cns)
   {
      if (cns.type() == CopyNumbersSection::types::Average)
      {
         cna_ = std::make_unique<CopyNumbersAvg>(world_, rules_, cns.interval());
         simulator_->add_extrnal_event(0, cna_.get());
         simulator_->add_reaction_recorder(cna_.get());
         if (!cns.file().empty())
         {
            if (!cfile_.is_open()) cfile_.open(cns.file(), std::fstream::out | std::fstream::trunc);
            cna_->set_output(cfile_);
         }
      }
      else
      {
         cni_ = std::make_unique<CopyNumbersInst>(world_, cns.interval());
         simulator_->add_extrnal_event(0, cni_.get());
         if (!cns.file().empty())
         {
            if (!cfile_.is_open()) cfile_.open(cns.file(), std::fstream::out | std::fstream::trunc);
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
         if (!pfile_.is_open()) pfile_.open(pps.file(), std::fstream::out | std::fstream::trunc);
         pp_->set_output(pfile_);
      }
   }

   void set_reactionrecord_section(const ReactionRecordSection &rrs)
   {
      rrec_ = std::make_unique<reaction_recorder_log>(rules_, model_);
      simulator_->add_reaction_recorder(rrec_.get());
      if (!rrs.file().empty())
      {
         if (!rfile_.is_open()) rfile_.open(rrs.file(), std::fstream::out | std::fstream::trunc);
         rrec_->set_output(rfile_);
      }
   }

   void set_reactionrule_section(bool first)
   {
      for (auto& section : settings_.getReactionRuleSections())
      {
         if (test_mode_setup(first, section.mode()))
            section.create_reaction_rule(model_, rules_);
      }
   }

   void add_particles_section(bool first)
   {
      for (auto& section : settings_.getParticlesSections())
      {
         if (test_mode_setup(first, section.mode()))
            section.add_particles_to_world(model_, world_, rng_);
      }
   }

   void set_progress_section(bool first, uint width)
   {
      if (first && (end_time_ > 0 || prep_time_ > 0))
      {
         progress_ = std::make_unique<Progress>(Progress(prep_time_ > 0 ? prep_time_ : end_time_, width));
         simulator_->add_extrnal_event(0, progress_.get());
      }
      // Need to redo the progressbar (or other external events) since simulator is reset after pre-simulation phase.
      if (!first && end_time_ > 0)
      {
         progress_ = std::make_unique<Progress>(Progress(end_time_, width));
         simulator_->add_extrnal_event(0, progress_.get());
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string file_;
   SimulatorSettings settings_;

   std::fstream cfile_, pfile_, rfile_;
   std::unique_ptr<CopyNumbersInst> cni_;
   std::unique_ptr<CopyNumbersAvg> cna_;
   std::unique_ptr<ParticlePositions> pp_;
   std::unique_ptr<reaction_recorder_log> rrec_;
   std::unique_ptr<Progress> progress_;          // progress bar indicator

};

// --------------------------------------------------------------------------------------------------------------------------------
