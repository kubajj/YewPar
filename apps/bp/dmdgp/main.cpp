/****************************************************************************************************
  Name:       MD-jeep
              the Branch & Prune algorithm for discretizable Distance Geometry - main program
  Author:     A. Mucherino, D.S. Goncalves, C. Lavor, L. Liberti, J-H. Lin, N. Maculan
  Sources:    ansi C
  License:    GNU General Public License v.3
  History:    May 01 2010  v.0.1    first release
              May 10 2014  v.0.2    method based on change of basis added for coordinate computation
                                    more robust check on input instances added
              Jul 28 2019  v.0.3.0  main adapted for BP with interval distances
                                    more efficient organization of the data structures
              Mar 21 2020  v.0.3.1  using new functions of "vertex" to verify the instance properties
                                    precomputing all triplets of reference vertices
              May 19 2020  v.0.3.2  introduction of MDfiles, possibility to select the method to run
              Apr 13 2022  v.0.3.2  patch

  Code adapted to work with YewPar by Jakub Jelinek and Blair Archibald
*****************************************************************************************************/
#include "bp.hpp"

double INFTY = 1.e+30;

int hpx_main(hpx::program_options::variables_map &opts)
{
  /*
  if (!opts.count("input-file")) {
    std::cerr << "You must provide an DIMACS input file with \n";
    hpx::finalize();
    return EXIT_FAILURE;
  }
  */
  FILE *input;
  OPTION op;
  INFORMATION info;
  size_t nlines, linelen, wordlen;
  char *line;
  int n0;

  // Check if the help option was provided
  // if (vm.count("help"))
  // {
  //   hpx::cout << desc_commandline << std::endl
  //             << "Note: When using -1, options -p and -P have the same effect" << std::endl;
  //   hpx::finalize();
  //   return EXIT_FAILURE;
  // }

  // hpx::program_options::notify(opts);

  auto inputFile = opts["mdfile"].as<std::string>();

  ParsedData data = parseFile(inputFile);
  std::vector<DataRecord> instance = readDataFile(data.file);

  auto start_time = std::chrono::steady_clock::now();

  auto skeletonType = opts["skeleton"].as<std::string>();

  auto overall_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time);

  hpx::cout << "cpu = " << overall_time.count() << std::endl;

  return hpx::finalize();
}

int main(int argc, char *argv[])
{
  hpx::program_options::options_description
      desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

  // clang-format off
  desc_commandline.add_options()
      ("help,h", "Print this help message")
      ( "skeleton-type",
        hpx::program_options::value<std::string>()->default_value("seq"),
        "Which skeleton to use: seq, depthbound, stacksteal, budget, or ordered"
      )
      ( "mdfile",
        hpx::program_options::value<std::string>()->required(),
        "Whhere to find the MDFile which contains specification of the instance and necessary fields"
      )
      ( "1", "The specified method stops at the first solution")
      ( "-l",  hpx::program_options::value<unsigned>(),
        "Specifies after how many solutions the method should stop"
      )
      ( "p", "Prints the best found solution in a text file")
      ( "P", "Prints all found solutions")
      ( "solutions,l",
        hpx::program_options::value<unsigned>(),
        "Specifies after how many solutions the method should stop"
      )
      ( "format, f",
        hpx::program_options::value<std::string>()->default_value("xyz"),
        "Specifies the output format (default is \"xyz\", may be changed to \"pdb\")"
      )
      ( "consec", "Verifies whether the consecutivity assumption is satisfied")
      ( "nomonitor", "Does not show the current layer number during the execution to improve performance")
      ("sym", hpx::program_options::value<unsigned>()->default_value(2), "Only one symmetric half of the tree is explored (argument may be 1 or 2)");
  // clang-format on

  hpx::register_startup_function(&Workstealing::Policies::SearchManagerPerf::registerPerformanceCounters);

  hpx::init_params args;
  args.desc_cmdline = desc_commandline;
  return hpx::init(argc, argv, args);
}
