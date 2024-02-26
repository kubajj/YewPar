/*

Simple application adapting the Bron-Kerbosh algoprithm to enumerate Maximal Cliques

*/

#include <iostream>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <set>
#include <chrono>
#include <memory>
#include <typeinfo>

#include <hpx/hpx_init.hpp>
#include <hpx/iostream.hpp>

#include <boost/serialization/access.hpp>

#include "DimacsParser.hpp"

#include "YewPar.hpp"

#include "skeletons/Seq.hpp"
#include "skeletons/DepthBounded.hpp"
#include "skeletons/StackStealing.hpp"
#include "skeletons/Ordered.hpp"
#include "skeletons/Budget.hpp"

#include "util/func.hpp"
#include "util/NodeGenerator.hpp"

struct BKNode
{
    friend class boost::serialization::access;
    std::set<int> R, P, X;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & R;
        ar & P;
        ar & X;
    }
};

struct SearchSpace
{
    std::map<int, std::set<int>> graph;

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & graph;
    }
};

struct GenNode : YewPar::NodeGenerator<BKNode, SearchSpace>
{
    std::set<int> R, P, X;
    int i, v;
    std::set<int>::iterator iter;
    std::reference_wrapper<const SearchSpace> space;

    // constructor
    GenNode(const SearchSpace &space, const BKNode &node) : space(std::cref(space))
    {
        // Body
        if (node.P.empty() && node.X.empty())
        {
            // R is maximal clique
            std::cout << "Found maximal clique" << std::endl;
            numChildren = 0;
        }
        else
        {
            R = node.R;
            P = node.P;
            X = node.X;
            numChildren = P.size();
            i = 0;
            iter = P.begin();
        }
    }

    // Return the next BKNode to look into
    BKNode next() override
    {
        if (iter == P.end())
        {
            std::cerr << "Wrong number of children or wrong iterator" << std::endl;
        }
        v = *iter;
        BKNode newNode;
        std::set<int> vSet = {v};
        auto neighbourSet = space.get().graph.find(v)->second;
        std::set_union(R.begin(), R.end(), vSet.begin(), vSet.end(), std::inserter(newNode.R, newNode.R.begin()));
        std::set_intersection(P.begin(), P.end(), neighbourSet.begin(), neighbourSet.end(), std::inserter(newNode.P, newNode.P.begin()));
        std::set_intersection(X.begin(), X.end(), neighbourSet.begin(), neighbourSet.end(), std::inserter(newNode.X, newNode.X.begin()));
        iter = P.erase(iter);
        X.insert(v);
        return newNode;
    }
};

struct CountSols : YewPar::Enumerator<Node, std::uint64_t>
{
    std::uint64_t count;
    CountSols() : count(0){};

    void accumulate(const BKNode &n) override
    {
        if (n.P.empty() && n.X.empty())
        {
            count++;
        }
    }

    void combine(const std::uint64_t &other) override
    {
        count += other;
    }

    std::uint64_t get() override { return count; }
};

int hpx_main(hpx::program_options::variables_map &opts)
{
    auto inputFile = opts["input-file"].as<std::string>();

    auto gFile = dimacs::read_dimacs(inputFile);

    SearchSpace space;
    space.graph = gFile.second;
    BKNode root;
    for (unsigned i = 0; i < gFile.first; ++i)
    {
        root.P.insert(i);
    }

    auto start_time = std::chrono::steady_clock::now();

    std::uint64_t count;
    if (skeleton == "seq")
    {
        YewPar::Skeletons::API::Params<> searchParameters;
        count = YewPar::Skeletons::Seq<
            GenNode,
            YewPar::Skeletons::API::Enumeration,
            YewPar::Skeletons::API::Enumerator<CountSols>,
            YewPar::Skeletons::API::DepthLimited>::search(space, root, searchParameters);
    }
    else
    {
        hpx::cout << "Invalid skeleton type: " << skeleton << std::endl;
        return hpx::finalize();
    }

    auto overall_time = std::chrono::duration_cast<
        std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time);

    hpx::cout << "cpu = " << overall_time.count() << std::endl;

    return hpx::finalize();
}

int main(int argc, char *argv[])
{
    hpx::program_options::options_description
        desc_commandline("Usage: " HPX_APPLICATION_STRING " [options]");
    // clang-format off
    desc_commandline.add_options()
    ( "skeleton",
      hpx::program_options::value<std::string>()->default_value("seq"),
      "Which skeleton to use: seq, depthbound, stacksteal, budget, or ordered"
      )
    ( "spawn-depth,d",
      hpx::program_options::value<std::uint64_t>()->default_value(0),
      "Depth in the tree to spawn at"
      )
    ( "backtrack-budget,b",
      hpx::program_options::value<unsigned>()->default_value(50),
      "Number of backtracks before spawning work"
      )
    ( "input-file,f",
      hpx::program_options::value<std::string>()->required(),
      "DIMACS formatted input graph"
      )
    ("discrepancyOrder", "Use discrepancy order for the ordered skeleton")
    ("chunked", "Use chunking with stack stealing")
    ("poolType",
     hpx::program_options::value<std::string>()->default_value("depthpool"),
     "Pool type for depthbounded skeleton");
    // clang-format on
    YewPar::registerPerformanceCounters();

    hpx::init_params args;
    args.desc_cmdline = desc_commandline;
    return hpx::init(argc, argv, args);
}
