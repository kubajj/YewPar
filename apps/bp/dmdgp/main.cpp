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

#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
#include <map>
#include <chrono>
#include <memory>
#include <typeinfo>

#include <hpx/hpx_init.hpp>
#include <hpx/iostream.hpp>

#include <boost/serialization/access.hpp>

#include "YewPar.hpp"

#include "skeletons/Seq.hpp"
#include "skeletons/DepthBounded.hpp"
#include "skeletons/StackStealing.hpp"
#include "skeletons/Ordered.hpp"
#include "skeletons/Budget.hpp"

#include "bp.hpp"

#include "util/func.hpp"
#include "util/NodeGenerator.hpp"

int hpx_main(hpx::program_options::variables_map & opts) {
  /*
  if (!opts.count("input-file")) {
    std::cerr << "You must provide an DIMACS input file with \n";
    hpx::finalize();
    return EXIT_FAILURE;
  }
  */
  OPTION::op;
  INFORMATION::info;

  // Check if the help option was provided
  if (vm.count("help")) {
    std::cout << desc_commandline << std::endl << "Note: When using -1, options -p and -P have the same effect" << std::endl;
    return hpx::finalize();
  }

  //hpx::program_options::notify(opts);

  auto inputFile = opts["mdfile"].as<std::string>();

  auto gFile = readMDfile(inputFile, &op, &info);

  // Order the graph (keep a hold of the map)
  std::map<int, int> invMap;
  auto graph = orderGraphFromFile<NWORDS>(gFile, invMap);

  auto spawnDepth = opts["spawn-depth"].as<std::uint64_t>();
  auto decisionBound = opts["decisionBound"].as<int>();

  auto start_time = std::chrono::steady_clock::now();

  // Initialise Root Node
  MCSol mcsol;
  mcsol.members.reserve(graph.size());
  mcsol.colours = 0;

  BitSet<NWORDS> cands;
  cands.resize(graph.size());
  cands.set_all();
  MCNode root = { mcsol, 0, cands };

  auto sol = root;
  auto skeletonType = opts["skeleton"].as<std::string>();
  if (skeletonType == "seq") {
    if (decisionBound != 0) {
      YewPar::Skeletons::API::Params<int> searchParameters;
      searchParameters.expectedObjective = decisionBound;

      sol = YewPar::Skeletons::Seq<GenNode,
                                   YewPar::Skeletons::API::Decision,
                                   YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                   YewPar::Skeletons::API::PruneLevel>
            ::search(graph, root, searchParameters);
    } else {
    sol = YewPar::Skeletons::Seq<GenNode,
                                 YewPar::Skeletons::API::Optimisation,
                                 YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                 YewPar::Skeletons::API::PruneLevel>
          ::search(graph, root);
    }
  } else if (skeletonType == "depthbounded") {
    if (decisionBound != 0) {
      YewPar::Skeletons::API::Params<int> searchParameters;
      searchParameters.expectedObjective = decisionBound;
      searchParameters.spawnDepth = spawnDepth;
      sol = YewPar::Skeletons::DepthBounded<GenNode,
                                           YewPar::Skeletons::API::Decision,
                                           YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                           YewPar::Skeletons::API::PruneLevel>
            ::search(graph, root, searchParameters);
    } else {
      YewPar::Skeletons::API::Params<int> searchParameters;
      searchParameters.spawnDepth = spawnDepth;
      auto poolType = opts["poolType"].as<std::string>();
      if (poolType == "deque") {
        sol = YewPar::Skeletons::DepthBounded<GenNode,
                                             YewPar::Skeletons::API::Optimisation,
                                             YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                             YewPar::Skeletons::API::PruneLevel,
                                             YewPar::Skeletons::API::DepthBoundedPoolPolicy<
                                               Workstealing::Policies::Workpool> >
            ::search(graph, root, searchParameters);
      } else {
        sol = YewPar::Skeletons::DepthBounded<GenNode,
                                             YewPar::Skeletons::API::Optimisation,
                                             YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                             YewPar::Skeletons::API::PruneLevel,
                                             YewPar::Skeletons::API::DepthBoundedPoolPolicy<
                                               Workstealing::Policies::DepthPoolPolicy> >
            ::search(graph, root, searchParameters);
      }
    }
  } else if (skeletonType == "stacksteal") {
    if (decisionBound != 0) {
      YewPar::Skeletons::API::Params<int> searchParameters;
      searchParameters.expectedObjective = decisionBound;
      searchParameters.stealAll = static_cast<bool>(opts.count("chunked"));
      sol = YewPar::Skeletons::StackStealing<GenNode,
                                             YewPar::Skeletons::API::Decision,
                                             YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                             YewPar::Skeletons::API::PruneLevel>
          ::search(graph, root, searchParameters);
    } else {
      YewPar::Skeletons::API::Params<int> searchParameters;
      searchParameters.stealAll = static_cast<bool>(opts.count("chunked"));
      sol = YewPar::Skeletons::StackStealing<GenNode,
                                             YewPar::Skeletons::API::Optimisation,
                                             YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                             YewPar::Skeletons::API::PruneLevel>
          ::search(graph, root, searchParameters);
    }
  } else if (skeletonType == "ordered") {
    YewPar::Skeletons::API::Params<int> searchParameters;
    searchParameters.spawnDepth = spawnDepth;
    if (opts.count("discrepancyOrder")) {
      sol = YewPar::Skeletons::Ordered<GenNode,
                                       YewPar::Skeletons::API::Optimisation,
                                       YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                       YewPar::Skeletons::API::DiscrepancySearch,
                                       YewPar::Skeletons::API::PruneLevel>
          ::search(graph, root, searchParameters);
    } else {
    sol = YewPar::Skeletons::Ordered<GenNode,
                                         YewPar::Skeletons::API::Optimisation,
                                         YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                         YewPar::Skeletons::API::PruneLevel>
          ::search(graph, root, searchParameters);
    }
  } else if (skeletonType == "budget") {
    if (decisionBound != 0) {
    YewPar::Skeletons::API::Params<int> searchParameters;
    searchParameters.backtrackBudget = opts["backtrack-budget"].as<unsigned>();
    searchParameters.expectedObjective = decisionBound;
    sol = YewPar::Skeletons::Budget<GenNode,
                                    YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                    YewPar::Skeletons::API::Decision,
                                    YewPar::Skeletons::API::PruneLevel>
        ::search(graph, root, searchParameters);
    } else {
      YewPar::Skeletons::API::Params<int> searchParameters;
      searchParameters.backtrackBudget = opts["backtrack-budget"].as<unsigned>();
      sol = YewPar::Skeletons::Budget<GenNode,
                                      YewPar::Skeletons::API::Optimisation,
                                      YewPar::Skeletons::API::BoundFunction<upperBound_func>,
                                      YewPar::Skeletons::API::PruneLevel>
          ::search(graph, root, searchParameters);
    }
  } else {
    hpx::cout << "Invalid skeleton type option. Should be: seq, depthbound, stacksteal or ordered" << std::endl;
    hpx::finalize();
    return EXIT_FAILURE;
  }

  auto overall_time = std::chrono::duration_cast<std::chrono::milliseconds>
    (std::chrono::steady_clock::now() - start_time);

  hpx::cout << "MaxClique Size = " << sol.size << std::endl;
  hpx::cout << "cpu = " << overall_time.count() << std::endl;

  return hpx::finalize();
}

int main(int argc, char* argv[]) {
  hpx::program_options::options_description
      desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");

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
      ("sym", hpx::program_options::value<unsigned>()->default_value(2), "Only one symmetric half of the tree is explored (argument may be 1 or 2)")
        );

  hpx::register_startup_function(&Workstealing::Policies::SearchManagerPerf::registerPerformanceCounters);

  hpx::init_params args;
  args.desc_cmdline = desc_commandline;
  return hpx::init(argc, argv, args);
}
