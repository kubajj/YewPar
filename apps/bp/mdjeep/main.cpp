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

double INFTY = 1.e+30;

int hpx_main(hpx::program_options::variables_map &opts)
{
    int i, n, n0, m, mexact;
    int fidx, verr;
    int it, flag;
    size_t nlines, linelen, wordlen;
    bool clique, smallsine;
    bool check_consec;
    double **X;
    double obj, cosine;
    char *errmsg;
    char *line;
    char *timestring;
    VERTEX *v;
    SEARCH S;
    OPTION op;
    INFORMATION info;
    unsigned long typelist;
    struct timeval t1, t2;
    FILE *input;

    // hpx::program_options::notify(opts);

    std::cout << "YewPar DDGP solver based on MDjeep" << std::endl;

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