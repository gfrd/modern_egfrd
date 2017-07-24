// --------------------------------------------------------------------------------------------------------------------------------

#include "RunGfrd.hpp"
#include <gsl/gsl_errno.h>
#include <csignal>
#include <exceptions.hpp>

// --------------------------------------------------------------------------------------------------------------------------------

void print_usage()
{
   std::cout << "====================================================================\n";
   std::cout << "=   enhanced Green Function Reaction Dynamics (eGFRD) Simulation   =\n";
   std::cout << "====================================================================\n\n";
   std::cout << "RunGfrd type [options]" << std::endl << std::endl;
   std::cout << "   type      Select simulation model, type can be:" << std::endl;
   std::cout << "      Equilibrium     * info here *" << std::endl;
   std::cout << "      PowerSpectrum   * info here *" << std::endl;
   std::cout << "      Reversible      * info here *" << std::endl;
   std::cout << "      mapK            * info here *" << std::endl;
   std::cout << "      Custom          : Load simulation mode from input file." << std::endl;
   std::cout << "      Resume          : Continue a previously saved simulator dump file." << std::endl << std::endl;
}

// --------------------------------------------------------------------------------------------------------------------------------

void gsl_error_handler(const char* reason, const char* file, int line, int error)
{
   std::string ex_msg = make_string() << "GSL_ERROR: " << gsl_strerror(error) << " in " << file << ":" << line << " - " << reason;
   throw  gfrd_exception(ex_msg);
}

// --------------------------------------------------------------------------------------------------------------------------------

volatile bool _abort = false;
void signal_handler(int signal)
{
   if (signal == SIGINT) _abort = true;
}

// --------------------------------------------------------------------------------------------------------------------------------

int main(int argc, char** argv)
{
   bool arg_err = false;
   std::unique_ptr<Simulation> sim(nullptr);

   gsl_set_error_handler(&gsl_error_handler);
   std::signal(SIGINT, signal_handler);

   getoptions args(argc, argv);
   for (size_t i = 0; i < args.size() && !arg_err; ++i)
   {
      if (i == 0)
      {
         if (args.option(i) == "Equilibrium") sim = std::make_unique<Equilibrium>();
         else if (args.option(i) == "PowerSpectrum") sim = std::make_unique<PowerSpectrum>();
         else if (args.option(i) == "Reversible") sim = std::make_unique<Reversible>();
         else if (args.option(i) == "Resume") sim = std::make_unique<Resume>();
         else if (args.option(i) == "Custom") sim = std::make_unique<Custom>();
         else if (args.option(i) == "mapK") sim = std::make_unique<MapK>();
         else arg_err = true;
      }
      else if (sim) arg_err |= sim->HandleCommandArguments(i, args);
   }

   if (arg_err || sim == nullptr)
   {
      print_usage();
      if (sim) sim->print_usage();
      else std::cout << "   [options] Set additional options when simulation type is specified." << std::endl << std::endl;
      return 1;
   }

   sim->set_abort(_abort);
   sim->Run();
   return sim->failed() ? 2 : _abort ? 3 : 0;
}

// --------------------------------------------------------------------------------------------------------------------------------