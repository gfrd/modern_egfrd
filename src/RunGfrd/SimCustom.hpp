#pragma once

// --------------------------------------------------------------------------------------------------------------------------------

class SimCustom : public Simulation
{
public:

   // --------------------------------------------------------------------------------------------------------------------------------

   explicit SimCustom() noexcept
   {
      world_size_ = 1e-7;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Custom"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   int HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "m" && args.isvalue_NP(i + 1)) maintenance_step_ = std::stoi(args.option(++i));
      else if (args.isparam(i) && args.option(i) == "mf" && args.isvalue_F(i + 1)) simstate_file_ = args.option(++i);
      else if (args.isparam(i) && args.option(i) == "e" && args.isvalue_D(i + 1)) end_time_ = std::stod(args.option(++i));
      else return Simulation::HandleCommandArguments(i, args);
      return -1;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      std::cout << "  RunGfrd -c,--custom [ options ]" << std::endl;
      std::cout << "        [-h,-?,--help]        Print command line usage information" << std::endl;
      std::cout << "        [-m N]                Maintenance every N steps\n";
      std::cout << "        [-mf file]            Maintenance output file\n";
      std::cout << "        [-e time]             End simulation after model time in seconds\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   void PrintSettings() override
   {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      auto time_local = std::localtime(&time);
      std::cout << std::setw(14) << "time local = " << std::asctime(time_local);

      std::cout << "\n";
      
      Simulation::PrintSettings();

      std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      // Construct your simulation model here
      s1 = model_.add_species_type(SpeciesType("A", model_.get_def_structure_type_id(), 1e-12, 1e-9));
      s2 = model_.add_species_type(SpeciesType("B", model_.get_def_structure_type_id(), 1e-12, 1e-9));
      s3 = model_.add_species_type(SpeciesType("C", model_.get_def_structure_type_id(), 1e-12, 1e-9));

      // Create the world and simulator
      Simulation::SetupSimulation();

      // Add particles (three layers)
      world_.throwInParticles(s1, 24, rng_, false, Vector3(0, 0, 0), Vector3(world_size_, world_size_/3, world_size_) );
      world_.throwInParticles(s2, 24, rng_, false, Vector3(0, world_size_/3, 0), Vector3(world_size_, 2*world_size_/3, world_size_));
      world_.throwInParticles(s3, 24, rng_, false, Vector3(0, 2*world_size_/3, 0), Vector3(world_size_, world_size_, world_size_));

      // Add some rules
      rules_.add_reaction_rule(ReactionRule(s1, 24, std::vector < SpeciesTypeID > {s2}));
      rules_.add_reaction_rule(ReactionRule(s2, 24, std::vector < SpeciesTypeID > {s3}));
      rules_.add_reaction_rule(ReactionRule(s3, 24, std::vector < SpeciesTypeID > {s1}));

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   SpeciesTypeID s1,s2,s3;
};

// --------------------------------------------------------------------------------------------------------------------------------
