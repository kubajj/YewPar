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

  std::cout << "YewPar DDGP solver based on MDjeep" << std::endl;

  auto inputFile = opts["mdfile"].as<std::string>();
  const char *inputFile_c = inputFile.c_str();
  input = NULL;
  input = fopen(inputFile_c, "r");
  if (input == NULL)
  {
    fprintf(stderr, "mdjeep: error while opening MDfile '%s'\n", inputFile_c);
    return hpx::finalize();
  };

  // reading the MDfile and setting up the main options and infos
  errmsg = readMDfile(input, &op, &info);

  // any error occurred while reading MDfile?
  if (errmsg != NULL)
  {
    fprintf(stderr, "%s\n", errmsg);
    free(errmsg);
    return hpx::finalize();
  };
  fclose(input);

  // MDfile status
  fprintf(stderr, "mdjeep: MDfile read, instance name '%s'\n", info.name);
  fprintf(stderr, "mejeep: selected method is '");
  if (info.method == 0)
    fprintf(stderr, "bp");
  else
    fprintf(stderr, "spg");
  fprintf(stderr, "'\n");
  if (info.method == 0)
    fprintf(stderr, "mdjeep: tolerance epsilon = %g, resolution = %5.2lf, maxtime = %ds\n", op.eps, op.r, op.maxtime);
  if (info.refinement == 1)
    fprintf(stderr, "mdjeep: selected refinement method is 'spg'\n");
  else
    fprintf(stderr, "mdjeep: no refinement method selected\n");

  // setting up other default values for options and infos
  // (the ones not included in the MDfile)
  op.print = 0;
  op.format = 0;
  op.allone = 0;
  op.symmetry = 0;
  op.monitor = true;
  op.be = 0.10;
  info.exact = false;
  info.consec = false;
  info.ncalls = 0;
  info.nspg = 0;
  info.nspgok = 0;
  info.nsols = 0;
  info.maxsols = 10;
  info.pruning = 0;
  info.best_sol = 0;
  info.best_mde = INFTY;
  info.best_lde = INFTY;
  check_consec = false;

  // Not checking the other input arguments
  // additional information is printed on the screen (other mdjeep options)
  if (op.print == 1)
    fprintf(stderr, "mdjeep: the best solution ");
  if (op.print == 2)
    fprintf(stderr, "mdjeep: all solutions ");
  if (op.print != 0)
  {
    fprintf(stderr, "will be printed in ");
    if (op.format == 0)
      fprintf(stderr, "XYZ");
    else
      fprintf(stderr, "PDB");
    fprintf(stderr, " format\n");
  };
  if (op.allone == 1)
    fprintf(stderr, "mdjeep: only one solution is requested by the user\n");
  if (info.maxsols != 10)
    fprintf(stderr, "mdjeep: limit on maximum number of solutions is set to %d\n", info.maxsols);
  if (op.symmetry != 0)
    fprintf(stderr, "mdjeep: only one symmetric half of the tree is explored: ");
  if (op.symmetry == 1)
    fprintf(stderr, "left-hand subtree\n");
  if (op.symmetry == 2)
    fprintf(stderr, "right-hand subtree\n");

  // opening input file
  input = fopen(info.filename, "r");
  if (input == NULL)
  {
    fprintf(stderr, "mdjeep: cannot open instance file '%s'\n", info.filename);
    return hpx::finalize();
  };

  // verifying the length of words and lines in the text file (for proper memory allocations)
  nlines = textFileAnalysis(input, info.sep, &wordlen, &linelen);
  if (nlines == 0 || wordlen == 0 || linelen == 0)
  {
    fprintf(stderr, "mdjeep: error while reading instance file: the file seems to be empty\n");
    return hpx::finalize();
  };

  // memory allocation for the array of chars containing the text file lines
  line = (char *)calloc(linelen + 1, sizeof(char));

  // verifying whether the instance file is valid (all lines contain the same list of data types)
  if (!isDistanceFileValid(input, info.sep, &typelist, linelen, line))
  {
    fprintf(stderr, "mdjeep: error while reading instance file: different lines contain different lists of data types\n");
    free(line);
    return hpx::finalize();
  };

  // counting the number of vertices in the instance file
  n = numberOfVerticesInFile(input, info.sep, info.format, &n0, linelen, line);
  if (n == 0)
  {
    fprintf(stderr, "mdjeep: error while reading instance file: it looks like the instance file does not respect the specified format\n");
    free(line);
    return hpx::finalize();
  };

  // memory allocation for the array of vertex structures
  v = (VERTEX *)calloc(n, sizeof(VERTEX));

  // loading the instance file in memory
  verr = readDistanceFile(input, info.sep, n, n0, v, info.format, linelen, line);
  fclose(input);
  free(line);
  if (verr != -1)
  {
    fprintf(stderr, "mdjeep: error while reading instance file: ");
    if (verr == -2)
      fprintf(stderr, "it looks like the file does not respect the specified format\n");
    else if (verr == -3)
      fprintf(stderr, "the presence of a distance from a vertex to itself was detected\n");
    else if (verr == -4)
      fprintf(stderr, "some vertex ranks in the interval [%d,%d] are missing\n", n0, n0 + n);
    else if (verr == -5)
      fprintf(stderr, "some lower bounds are strictly greater than the corresponding upper bounds\n");
    else
      fprintf(stderr, "vertex with rank %d was found for the second time but with a different set of attributes\n", n0 + verr);
    free(v);
    return hpx::finalize();
  };

  // counting the number of distances
  m = totalNumberOfDistances(n, v);
  if (info.method == 0 && m < 3 * (n - 2))
  {
    fprintf(stderr, "mdjeep: error: not enough distances to perform discretization necessary to execute bp method\n");
    free(v);
    return hpx::finalize();
  };

  // counting the number of exact distances
  mexact = totalNumberOfExactDistances(n, v, op.eps);
  if (info.method == 0 && mexact < 2 * (n - 3) + 3)
  {
    fprintf(stderr, "mdjeep: error: not enough exact distances to perform discretization necessary to execute bp method\n");
    fprintf(stderr, "               a distance [lb,ub] is considered as exact if ub-lb < tolerance eps\n");
    free(v);
    return hpx::finalize();
  };

  // printing instance details
  fprintf(stderr, "mdjeep: instance file '%s' read: %d vertices / %d distances\n", info.filename, n, m);
  /*
    Main body
  */

  auto start_time = std::chrono::steady_clock::now();

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