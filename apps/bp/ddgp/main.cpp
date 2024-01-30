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

struct DDGPNode
{
    friend class boost::serialization::access;
    int id;
    int n_vertices;
    std::array<double, 12> qi;
    DDGPSol sol;

    int getObj() const
    {
        return id;
    }

    template <class Archive>
    void serialize(Archive &ar, const unsigned int version)
    {
        ar & id;
        ar & qi;
        ar & sol;
    }
};

struct CountSols : YewPar::Enumerator<DDGPNode, std::uint64_t>
{
    std::uint64_t count;
    CountSols() : count(0){};

    void accumulate(const DDGPNode &node) override
    {
        if (node.id == node.n_vertices)
            count++;
    }

    void combine(const std::uint64_t &other) override
    {
        count += other;
    }

    std::uint64_t get() override { return count; }
};

struct GenNode : YewPar::NodeGenerator<DDGPNode, DDGPMaps>
{

    // constructor
    GenNode(const DDGPMaps &maps, const DDGPNode &node)
    {
        // Body
        // Calculate Qi' - qi1 and Qi'' - qi2
        numChildren = 2;
    }

    // Return the next DDGPNode to look into
    DDGPNode next() override
    {
        DDGPNode nextNode;
        return nextNode;
    }
};

int hpx_main(hpx::program_options::variables_map &opts)
{

    // hpx::program_options::notify(opts);

    std::cout << "YewPar DDGP solver" << std::endl;

    auto inputFile = opts["mdfile"].as<std::string>();

    int n_vertices;

    ParsedData data = parseFile(inputFile);

    std::cout << "MDFile read, instance name " << data.file << std::endl;

    // Note(kubajj): n_vertices is 1..n (not 0 based)
    std::vector<DataRecord> instance = readDataFile(data.file, n_vertices);
    VertexReferences refs;
    DistanceMap distanceMap = createDataRecordMap(instance, refs);

    hpx::cout << "Number of vertices in the file: " << n_vertices << std::endl;

    auto start_time = std::chrono::steady_clock::now();

    /*
      Main body
    */

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
      ( "skeleton",
        hpx::program_options::value<std::string>()->default_value("seq"),
        "Which skeleton to use: only seq possible for now"
        )
      ( "mdfile",
        hpx::program_options::value<std::string>()->required(),
        "Where to find the MDFile which contains specification of the instance and necessary fields"
      );
    // clang-format on

    YewPar::registerPerformanceCounters();

    hpx::init_params args;
    args.desc_cmdline = desc_commandline;
    return hpx::init(argc, argv, args);
}