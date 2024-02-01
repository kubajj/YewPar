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
  int n;
  VERTEX *v;
  double **X;
  SEARCH S;
  OPTION op;
  INFORMATION info;

  int i, n0, m, mexact;
  int fidx, verr;
  int it, flag;
  size_t nlines, linelen, wordlen;
  bool clique, smallsine;
  bool check_consec;
  double obj, cosine;
  char *errmsg;
  char *line;
  unsigned long typelist;
  struct timeval t1, t2;
  FILE *input;

  std::cout << "YewPar DDGP solver based on MDjeep" << std::endl;

  auto inputFile = opts["mdfile"].as<std::string>();

  Config config = mdjeep_main(inputFile);
  if (config.error)
    return hpx::finalize();
  n = config.n;
  v = config.v;
  X = config.X;
  S = config.S;
  op = config.op;
  info = *(config.info);

  /*
    Main body
  */

  auto start_time = std::chrono::steady_clock::now();

  // calling method bp
  if (info.method == 0)
  {
    fprintf(stderr, "mdjeep: bp is exploring the search tree ... ");
    if (op.monitor)
    {
      fprintf(stderr, "layer ");
      for (i = 0; i < info.ndigits; i++)
        fprintf(stderr, " ");
    };
    gettimeofday(&t1, 0);
    if (info.exact)
      bp_exact(0, n, v, X, S, op, &info);
    else
      bp(0, n, v, X, S, op, &info);
    gettimeofday(&t2, 0);
    fprintf(stderr, "\n");
  };

  // printing the result found by bp
  if (info.method == 0)
  {
    if (t2.tv_sec - t1.tv_sec > op.maxtime)
      fprintf(stderr, "mdjeep: bp stopped because the maxtime was reached\n");
    fprintf(stderr, "mdjeep: %d solutions found by bp method", info.nsols);
    if (info.nsols == info.maxsols)
      fprintf(stderr, " (max %d)", info.maxsols);
    fprintf(stderr, "\n");
    fprintf(stderr, "mdjeep: %d branches were pruned\n", info.pruning);
    if (!info.exact)
      fprintf(stderr, "mdjeep: %d calls to spectral projected gradient (%d successful)\n", info.nspg, info.nspgok);
    if (info.nsols > 0)
      fprintf(stderr, "mdjeep: best solution #%d: LDE = %10.8lf, MDE = %10.8lf\n", info.best_sol, info.best_lde, info.best_mde);
  };

  auto overall_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time);

  hpx::cout << "cpu = " << overall_time.count() << std::endl;
  // freeing memory
  freeVector(S.memory);
  freeVector(S.Dy);
  freeVector(S.Yy);
  freeVector(S.Zy);
  freeMatrix(3, S.DX);
  freeMatrix(3, S.YX);
  freeMatrix(3, S.ZX);
  freeMatrix(3, S.Xp);
  freeMatrix(3, S.gXp);
  freeMatrix(3, S.gX);
  freeMatrix(3, S.sX);
  freeVector(S.yp);
  freeVector(S.gyp);
  freeVector(S.y);
  freeVector(S.gy);
  freeVector(S.sy);
  freeMatrix(3, S.pX);
  freeMatrix(3, S.lX);
  freeMatrix(3, S.uX);
  free(S.sym);
  freeMatrix(3, X);
  free(info.name);
  free(info.filename);
  if (info.method == 0)
    free(S.refs);
  if (op.print != 0)
    free(info.output);
  freeVertex(n, v);
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