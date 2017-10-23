// --------------------------------------------------------------------------------------------------------------------------------

#include <gsl/gsl_errno.h>
#include <csignal>
#include <exceptions.hpp>
#include "getoptions.hpp"
#include "Simulation.hpp"
#include "SimResume.hpp"
#include "SimModel.hpp"
#include "SimCustom.hpp"
#include "gfrd_compat.hpp"

// --------------------------------------------------------------------------------------------------------------------------------

void print_header()
{
   std::cout << "eGFRD (modern), NWO-I AMOLF, 2017, Amsterdam, The Netherlands, www.grfd.org" << std::endl;
}

void print_version()
{
   print_header();
   std::cout << " version: " << GFRD_VERSION_MAJOR << "." << GFRD_VERSION_MINOR << std::endl;
   std::cout << " build: "<< GFRD_VERSION_BUILD << std::endl;
   std::cout << " platform: " << (sizeof(std::nullptr_t)==8 ? "x64" : "x86");
#if _DEBUG
   std::cout << " DEBUG";
#endif
   std::cout << std::endl << " matrixspace: " << CompileConfigSimulator::MatrixCellsX << "x" << CompileConfigSimulator::MatrixCellsY << "x" << CompileConfigSimulator::MatrixCellsZ << std::endl;
   std::cout << std::endl;
}

void print_usage()
{
   std::cout << "  1) RunGfrd <path-to-model> [ options ]" << std::endl;
   std::cout << "  2) RunGfrd -r,--resume <path-to-simstate> [ options ]" << std::endl;
   std::cout << "  3) RunGfrd -c,--custom [ options ]" << std::endl;
   std::cout << "        [-h,-?,--help]        Print command line usage information" << std::endl;
   std::cout << "        [-v,--version]        Print version/build information" << std::endl  << std::endl;
   std::cout << "1) Start simulation described in model-file." << std::endl;
   std::cout << "2) Resume a simulation from state-file." << std::endl;
   std::cout << "3) Start simulation based on user-code in class SimCustom" << std::endl;
   std::cout << std::endl << std::endl;
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
   std::unique_ptr<Simulation> sim(nullptr);

   try 
   {
      int arg_err = -1;
      getoptions args(argc, argv);
      for (size_t i = 0; i < args.size() && arg_err==-1; ++i)
      {
         if (sim == nullptr)
         {
            if (!args.isparam(i))
            {
               sim = std::make_unique<SimModel>(args.option(i));        // path-to-model
               continue;
            }

            if (args.option(i) == "r" || args.option(i) == "-resume")
            { 
               if (args.isvalue_F(i + 1))
               {
                  sim = std::make_unique<SimResume>(args.option(++i));     // path-to-sim-state
                  continue; 
               }
            }

            if (args.option(i) == "c" || args.option(i) == "-custom")
            {
               sim = std::make_unique<SimCustom>();            // user-code simulation
               continue;
            }

        }

        if (args.option(i) == "h" || args.option(i) == "?" || args.option(i) == "-help")
        {
           std::cout << "Usage:" << std::endl << std::endl;
           if (sim) sim->print_usage(); else print_usage();
           return 1;
        }

        if (args.option(i) == "v" || args.option(i) == "-version")
        {
            print_version(); return 1;
        }

        if (sim) arg_err = sim->HandleCommandArguments(i, args); else arg_err = static_cast<int>(i);
      }

      if (arg_err!=-1 || sim == nullptr)
      {
         print_header();
         if (arg_err!=-1) 
            std::cout << "ERROR: Unknown or invalid argument: " << (args.isparam(arg_err) ? "-" : "") << args.option(arg_err) << std::endl;
         else 
            std::cout << "ERROR: Please specify your intentions." << std::endl;
         std::cout << "use --help argument to print usage information." << std::endl;
         return 1;
      }
   }
   catch (std::runtime_error ex)
   {
      Log("RunGfrd").fatal() << ex.what();
      return 2;
   }

   gsl_set_error_handler(&gsl_error_handler);
   std::signal(SIGINT, signal_handler);
   sim->set_abort(_abort);

   sim->Run();
   return sim->failed() ? 2 : _abort ? 3 : 0;
}

// --------------------------------------------------------------------------------------------------------------------------------