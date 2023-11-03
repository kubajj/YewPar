#include <hpx/hpx_init.hpp>
#include <hpx/iostream.hpp>

#include <boost/serialization/access.hpp>

#include "YewPar.hpp"

#include "skeletons/Seq.hpp"
#include "skeletons/DepthBounded.hpp"
#include "skeletons/StackStealing.hpp"
#include "skeletons/Ordered.hpp"
#include "skeletons/Budget.hpp"

#include "util/func.hpp"
#include "util/NodeGenerator.hpp"

#include "bp.hpp"

int hpx_main(hpx::program_options::variables_map &opts)
{
  /*
  if (!opts.count("input-file")) {
    std::cerr << "You must provide an DIMACS input file with \n";
    hpx::finalize();
    return EXIT_FAILURE;
  }
  */

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

  int n_vertices;

  ParsedData data = parseFile(inputFile);
  std::vector<DataRecord> instance = readDataFile(data.file, n_vertices);
  std::map<std::pair<int, int>, DataRecord *> references = createDataRecordMap(instance);

  // Create empty maps for theta and omega
  std::map<std::pair<int, int>, double> thetaMap;
  std::map<std::pair<int, int>, double> omegaMap;

  // Calculate the angles θ among consecutive triplets of vertices
  calculateAnglesForVertices(references, n_vertices, thetaMap, omegaMap);

  hpx::cout << "n: " << n_vertices << " len(thetaMap): " << thetaMap.size() << " len(omegaMap): " << omegaMap.size() << std::endl;

  auto start_time = std::chrono::steady_clock::now();

  auto skeletonType = opts["skeleton"].as<std::string>();

  /*
    Main body
  */
  if (skeletonType != "seq")
  {
    hpx::cout << "Invalid skeleton type option. Only seq implemented so far." << std::endl;
    hpx::finalize();
    return EXIT_FAILURE;
  }

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

  YewPar::registerPerformanceCounters();

  hpx::init_params args;
  args.desc_cmdline = desc_commandline;
  return hpx::init(argc, argv, args);
}
