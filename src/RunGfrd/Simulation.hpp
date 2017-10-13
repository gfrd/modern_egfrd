#ifndef SIMULATION_HPP
#define SIMULATION_HPP

// --------------------------------------------------------------------------------------------------------------------------------

#include <chrono>
#include <EGFRDSimulator.hpp>
#include <Progress.hpp>
#include <Persistence.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

class Simulation
{
public:

   // --------------------------------------------------------------------------------------------------------------------------------
   //
   // Order of method calls:
   //
   //  - HandleCommandArguments(..);                                 one call for each argument on the command line.
   //  - if arguments contain error: print_usage() and terminate
   //
   //  - Run();
   //    - PrepareSimulation();                                     calculate parameters and runtime settings
   //    - PrintSettings();                                         display
   //    - SetupSimulation();                                       create model/world/simulator
   //      - if (prep_time>0)
   //         - LoopSimulatorUntil(preptime)                        let system reach equilibrium (if required)
   //         - reset simulator
   //         - PostPreSimulation();                                setup true experiment, e.q. add reaction rules
   //    - LoopSimulatorUntil(end_time)                             
   //    - PostSimulation();                                        post processing and show results
   //
   // From within the inner simulation loops, SimStep() is called (processing per step)
   // and the maintenance() is called every maintenance_step iterations.
   // If maintenance fails the simulation is stopped.
   //
   // --------------------------------------------------------------------------------------------------------------------------------

   explicit Simulation() : world_size_(1e-7), prep_time_(0.0), end_time_(0.0), seed_(0), show_progress_(false), maintenance_step_(10000), simstate_file_(), failed_(false), abort_(nullptr) { }

   virtual ~Simulation() { }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual std::string name() const { return "Simulation"; }
   double end_time() const { return end_time_; }
   double prepare_time() const { return prep_time_; }
   double time() const { return simulator_ ? simulator_->time() : 0.0; }
   size_t num_steps() const { return simulator_ ? simulator_->num_steps() : 0; }
   void set_abort(volatile bool& abort) { abort_ = &abort; }
   bool failed() const { return failed_; }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual bool HandleCommandArguments(size_t& i, const getoptions& args)
   {
      if (!args.isparam(i)) return true; // expect a parameter
      if (args.option(i) == "ws" && args.isvalue(i + 1)) world_size_ = args.value(++i);
      else if (args.option(i) == "e" && args.isvalue(i + 1)) end_time_ = args.value(++i);
      else if (args.option(i) == "p" && args.isvalue(i + 1)) prep_time_ = args.value(++i);
      else if (args.option(i) == "s" && args.isvalue(i + 1)) seed_ = std::stoul(args.option(++i), nullptr, 0); // auto detect, hex/oct/dec
      else if (args.option(i) == "m" && args.isvalue(i + 1)) maintenance_step_ = args.number(++i);
      else if (args.option(i) == "ms" && args.isvalue(i + 1)) simstate_file_ = args.option(++i);
      else if (args.option(i) == "sp0") show_progress_ = false;
      else if (args.option(i) == "sp1") show_progress_ = true;
      else return true;
      return false;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void print_usage()
   {
      std::cout << "      -ws x.xx         World Size (cube side length in meter)\n";
      std::cout << "      -p x.xx          Prepare Time (equalize simulation for x.xx seconds)\n";
      std::cout << "      -e x.xx          End Time (stop simulation after seconds)\n";
      std::cout << "      -s nnnn          Random number seed (auto detect decimal/octal/hex value)\n";
      std::cout << "      -m n             Run maintenance every N simulation steps\n";
      std::cout << "      -ms filename     Store simulator state at maintenance step\n";
      std::cout << "      -sp0/-sp1        Show Progress Bar (0=off [default], 1=on) DO specify an end_time!\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   // Starts the EGFRD simulation (split into setup-phase, preparation run and main simulation run).
   void Run()
   {
      Log("RunGfrd").info() << "Setup " << name() << " simulation";

      try
      {
         if (!PrepareSimulation()) return;
         PrintSettings();
         if (!SetupSimulation()) return;
      }
      catch (std::runtime_error ex)
      {
         Log("RunGfrd").fatal() << ex.what();
         failed_ = true;
         return;
      }

      auto begin = std::chrono::high_resolution_clock::now();
      if (prep_time_ > 0)
      {
         Log("RunGfrd").info() << "Start of pre-simulation";
         LoopSimulatorUntil(prep_time_, true);
         if (!failed_ && (abort_ ? !*abort_ : true))
         {
            auto now = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::duration<float>>(now - begin).count();
            Log("RunGfrd").info() << "Pre-simulation took  " << std::fixed << std::setprecision(3) << elapsed << " seconds / " << num_steps() << " steps.";
            simulator_->reset();
            PostPreSimulation();
         }
      }

      if (!failed_ && (abort_ ? !*abort_ : true))
      {
         begin = std::chrono::high_resolution_clock::now();
         Log("RunGfrd").info() << "Start of simulation";
         LoopSimulatorUntil(end_time_);
         if (!failed_) PostSimulation();
      }

      auto now = std::chrono::high_resolution_clock::now();
      auto elapsed = std::chrono::duration_cast<std::chrono::duration<float>>(now - begin).count();
      Log("RunGfrd").info() << "Simulation " << (failed_ ? "failed after" : (abort_ ? *abort_ : false) ? "aborted after" : "took") << " : " << std::fixed << std::setprecision(3) << elapsed << " seconds / " << num_steps() << " steps.";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

protected:

   void LoopSimulatorUntil(double time, bool prerun = false)
   {
      try
      {
         while (!failed_ && (abort_ ? !*abort_ : true) && (time > 0.0 ? simulator_->time() < time : true))
         {
            if (maintenance_step_ != 0 && num_steps() % maintenance_step_ == 0)
            {
               failed_ = !maintenance();
               if (failed_) break;
            }

            simulator_->step();
            SimStep(prerun);
         }
      }
      catch (std::runtime_error ex)
      {
         failed_ = true;
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual bool maintenance()
   {
      try
      {
         simulator_->check();
      }
      catch (gfrd_exception ex)
      {
         Log("RunGfrd").fatal() << "Check failed in simulation, at time: " << time() << ", step: " << num_steps() << ", fault: " << ex.what();
         simulator_->dump(make_string() << "sim_check_failed_atstep" << std::setfill('0') << std::setw(10) << num_steps() << ".log", false);
         return false;
      }
      if (!simstate_file_.empty())
      {
         try
         {
            Persistence p;
            p.store(simstate_file_);
            p.store_egfrd(*simulator_);
         }
         catch (std::runtime_error ex)
         {
            Log("RunGfrd").warn() << "Failed to store simulator state to file! fault: " << ex.what();
         }
      }
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void SimStep(bool prerun)
   {
      UNUSED(prerun);
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual bool PrepareSimulation()
   {
      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void PrintSettings()
   {
      auto now = std::chrono::system_clock::now();
      auto time = std::chrono::system_clock::to_time_t(now);
      auto time_local = std::localtime(&time);
      std::cout << std::setw(14) << "time local = " << std::asctime(time_local);
      std::cout << std::setw(14) << "world size = " << world_size_ << " [m]\n";
      auto ms = world_.matrix_size();
      std::cout << std::setw(14) << "volume = " << std::setprecision(3) << 1e3 * (world_size_ * world_size_ * ms[1] / ms[0] * world_size_ * ms[2] / ms[0]) << " [L]\n";
      std::cout << std::setw(14) << "matrix size = " << ms[0] << "x" << ms[1] << "x" << ms[2] << "\n";
      if (seed_) std::cout << std::setw(16) << "seed = 0x" << std::setw(8) << std::setfill('0') << std::hex << std::uppercase << seed_ << std::setfill(' ') << "\n" << std::dec;
      if (prep_time_ > 0.0) std::cout << std::setw(14) << "prep_time = " << prep_time_ << " [s]\n";
      if (end_time_ > 0.0) std::cout << std::setw(14) << "end_time = " << end_time_ << " [s]\n";
      if (maintenance_step_ > 0) std::cout << std::setw(14) << "maintenance = " << maintenance_step_ << "\n";
      if (!simstate_file_.empty()) std::cout << std::setw(14) << "simstatefile = " << simstate_file_ << "\n";
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual bool SetupSimulation()
   {
      if (seed_) rng_.seed(seed_);

      world_.initialize(world_size_, model_);

      simulator_ = std::make_unique<EGFRDSimulator>(world_, rules_, rng_);

      if (show_progress_ && (end_time_ > 0 || prep_time_ > 0))
      {
         progress_ = std::make_unique<Progress>(Progress(prep_time_ > 0 ? prep_time_ : end_time_, 80));
         simulator_->add_extrnal_event(0, progress_.get());
      }

      return true;
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void PostPreSimulation()
   {
      // Need to redo the progressbar (or other external events) since simulator is reset after pre-simulation phase.
      if (show_progress_ && end_time_ > 0)
      {
         progress_ = std::make_unique<Progress>(Progress(end_time_, 80));
         simulator_->add_extrnal_event(0, progress_.get());
      }
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   virtual void PostSimulation()
   {
   }

   // --------------------------------------------------------------------------------------------------------------------------------

   double world_size_;                           // The edge length of the cube shaped world [m]
   double prep_time_;                            // The simulation prepare time (zero for no preparation run);
   double end_time_;                             // The simulation end time
   int seed_;                                    // Set seed of the random number generator (if non-zero)
   bool show_progress_;                          // The progress indicator in console is enabled when set to true
   Model model_;                                 // The model contains the species and structures
   World world_;                                 // The world contains all structures and particles
   ReactionRuleCollection rules_;                // The reaction rules between the structures and particles
   RandomNumberGenerator rng_;                   // The random number generator
   std::unique_ptr<EGFRDSimulator> simulator_;   // The EGFRD simulator
   std::unique_ptr<Progress> progress_;          // progress bar indicator
   size_t maintenance_step_;                      // run internal maintenance cycle every n sim steps
   std::string simstate_file_;                    // filename to store simulator state (in maintenance step)
   bool failed_;                                  // simulation failed (check or exception)
   volatile bool* abort_;                         // Ctrl-C abort of simulation loop

   // --------------------------------------------------------------------------------------------------------------------------------

};

#endif /* SIMULATION_HPP */
