#ifndef SIMRESUME_HPP
#define SIMRESUME_HPP

// --------------------------------------------------------------------------------------------------------------------------------

class SimResume : public Simulation
{
public:

   // --------------------------------------------------------------------------------------------------------------------------------

   explicit SimResume(const std::string& filename) noexcept : file_(filename) {}

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Resume"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   int HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "m" && args.isvalue(i + 1)) maintenance_step_ = std::stoi(args.option(++i));
      else if (args.isparam(i) && args.option(i) == "e" && args.isvalue(i + 1)) end_time_ = std::stod(args.option(++i));
      else return Simulation::HandleCommandArguments(i, args);
      return -1;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      std::cout << "  RunGfrd -r/--resume <path-to-simstate> [ options ]" << std::endl;
      std::cout << "        [-h,-?,--help  ]      Print command line usage information" << std::endl;
      std::cout << "        [-m N]                Maintenance every N steps\n";
      std::cout << "        [-s file]             Maintenance output file\n";
      std::cout << "        [-e d]                Stop simulation at d.dd seconds simulation time\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   void PrintSettings() override
   {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      auto time_local = std::localtime(&time);
      std::cout << std::setw(14) << "time local = " << std::asctime(time_local);
      std::cout << std::setw(14) << "resume file = " << file_ << "\n";

      std::cout << "\n";
      
      Simulation::PrintSettings();

      std::cout << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool SetupSimulation() override
   {
      Simulation::SetupSimulation();
      Persistence::callback_fn cb(nullptr); // TODO get lamda in here that maps CustomActions to members

      Persistence p;
      if (!p.retreive(file_, cb)) { Log("RunGfrd").fatal() << "Failed to load resume file!"; return false; }
      p.retreive_egfrd(*simulator_);

      // check after load!
      if (!maintenance()) { Log("RunGfrd").fatal() << "Maintenance check failed!"; return false; }

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   std::string file_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /*RESUME_HPP*/