#ifndef RESUME_HPP
#define RESUME_HPP

// --------------------------------------------------------------------------------------------------------------------------------

class Resume : public Simulation
{
public:

   explicit Resume() noexcept : resumefile_() {}

   // --------------------------------------------------------------------------------------------------------------------------------

   std::string name() const override { return "Resume"; }

   // --------------------------------------------------------------------------------------------------------------------------------

   bool HandleCommandArguments(size_t& i, const getoptions& args) override
   {
      if (args.isparam(i) && args.option(i) == "f" && args.isvalue(i + 1)) resumefile_ = args.option(++i);
      else return true; // Simulation::HandleCommandArguments(i, args);
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   void print_usage() override
   {
      std::cout << "      -f filename      simState file to resume\n";
      //Simulation::print_usage();
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   void PrintSettings() override
   {
      Simulation::PrintSettings();
      std::cout << std::setw(14) << "input file = " << resumefile_ << "\n";
   }

   bool SetupSimulation() override
   {
      if (resumefile_.empty()) { Log("RunGfrd").fatal() << "No resume file! Specify with -f command option"; return false; }

      Simulation::SetupSimulation();
      Persistence::callback_fn cb(nullptr); // TODO get lamda in here that maps CustomActions to members

      Persistence p;
      if (!p.retreive(resumefile_, cb)) { Log("RunGfrd").fatal() << "Failed to load resume file!"; return false; }
      p.retreive_egfrd(*simulator_);

      // check after load!
      if (!maintenance()) { Log("RunGfrd").fatal() << "Maintenance check failed!"; return false; }

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

private:
   std::string resumefile_;
};

// --------------------------------------------------------------------------------------------------------------------------------

#endif /*RESUME_HPP*/