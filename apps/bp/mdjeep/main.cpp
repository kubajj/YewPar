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

struct SearchNode
{
  friend class boost::serialization::access;
  int i;
  int n_vertices;
  double **X;
  SEARCH S;

  int getObj() const
  {
    return i;
  }

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version)
  {
    ar & i;
    ar & n_vertices;
    ar & X;
    ar & S;
  }
};

struct CountSols : YewPar::Enumerator<SearchNode, std::uint64_t>
{
  std::uint64_t count;
  CountSols() : count(0){};

  void accumulate(const SearchNode &node) override
  {
    if (node.i == node.n_vertices)
      count++;
  }

  void combine(const std::uint64_t &other) override
  {
    count += other;
  }

  std::uint64_t get() override { return count; }
};

struct GenNode : YewPar::NodeGenerator<SearchNode, SearchSpace>
{
  double **X1, **X2;
  SEARCH S1, S2;
  bool firstIsPruned;
  bool first;
  int i, n_vertices;

  // constructor
  GenNode(const SearchSpace &space, const SearchNode &node)
  {
    i = node.i;
    n_vertices = node.n_vertices;
    if (i == node.n_vertices)
    {
      numChildren = 0;
    }
    else
    {
      bool pruned;
      int j, k;
      int it, nb;
      double U[9];
      double cdist;
      double cTheta, sTheta;
      double cosOmega00, cosOmega01;
      double sinOmega00, sinOmega01;
      double lomega0, uomega0;
      double lomega1, uomega1;
      Omega *current;
      omegaList omegaL;
      REFERENCE *r1, *r2, *r3;
      // i = node.id;
      // updating BP call counter
      space.info->ncalls++;
      it = 0;

      // reference vertices
      r3 = node.S.refs[i].r3;
      r2 = node.S.refs[i].r2;
      r1 = node.S.refs[i].r1;
      cdist = lowerBound(r1);

      // theta angle ("bond" angles)
      cTheta = costheta(otherVertexId(r2), otherVertexId(r1), i, space.v, node.X);
      sTheta = sqrt(1.0 - cTheta * cTheta);

      // generating U matrix (only once)
      UMatrix(otherVertexId(r3), otherVertexId(r2), otherVertexId(r1), i, node.X, U);

      // omega angle (torsion angles)
      nb = 2;
      cosOmega00 = cosomega(otherVertexId(r3), otherVertexId(r2), otherVertexId(r1), i, space.v, node.X, 0.0, space.op.eps);
      cosOmega01 = cosomega(otherVertexId(r3), otherVertexId(r2), otherVertexId(r1), i, space.v, node.X, 1.0, space.op.eps);
      if (cosOmega00 == -2.0 || cosOmega01 == -2.0)
        return; // infeasibility already detected
      sinOmega00 = sqrt(1.0 - cosOmega00 * cosOmega00);
      sinOmega01 = sqrt(1.0 - cosOmega01 * cosOmega01);
      lomega0 = atan2(+sinOmega00, cosOmega00);
      uomega0 = atan2(+sinOmega01, cosOmega01);
      lomega1 = atan2(-sinOmega00, cosOmega00);
      uomega1 = atan2(-sinOmega01, cosOmega01);
      // if the two omega intervals are adjacent, we can consider the union
      if (i > 3)
      {
        if (fabs(uomega0 - lomega1) < space.op.eps)
        {
          nb = 1;
          uomega0 = uomega1;
        }
        else if (fabs(uomega1 - lomega0) < space.op.eps)
        {
          nb = 1;
          lomega0 = lomega1;
        };
      };
      // if the layer is symmetric, it is not necessary to consider the entire omega intervals
      if (node.S.sym[i])
      {
        lomega0 = 0.5 * (lomega0 + uomega0);
        uomega0 = lomega0;
        if (nb == 2)
        {
          lomega1 = 0.5 * (lomega1 + uomega1);
          uomega1 = lomega1;
        };
      };

      // initializing omega list
      omegaL = initOmegaList(lomega0, uomega0);
      if (nb == 2)
        attachNewOmegaInterval(firstOmegaInterval(omegaL), lomega1, uomega1);

      // verifying the "arclength" of every arc wrt the given resolution parameter
      splitOmegaIntervals(firstOmegaInterval(omegaL), cdist, space.op.r);

      // counting total number of omega intervals (necessary only at layer 3)
      if (i == 3)
        nb = numberOfOmegaIntervals(firstOmegaInterval(omegaL));

      // starting point for iterating over omega angles (it depends on op.symmetry)
      if (space.op.symmetry < 2)
        current = firstOmegaInterval(omegaL);
      else
        current = lastOmegaInterval(omegaL);

      X1 = copy(node.X, 3, n_vertices);
      X2 = copy(node.X, 3, n_vertices);
      S1 = copySearch(node.S, n_vertices, space.m);
      S2 = copySearch(node.S, n_vertices, space.m);

      numChildren = 0;
      pruned = prepare_branch(i, current, it, nb, cdist, cTheta, sTheta, U, r1, r3, X1, space.v, S1, space.op, space.info);
      if (!pruned)
      {
        // Assign X1
        numChildren++;
        firstIsPruned = false;
      }
      else
      {
        firstIsPruned = true;
      }
      if (omegaIntervalHasNextAlongDirection(current, space.op.symmetry < 2))
      {
        current = omegaIntervalNextAlongDirection(current, space.op.symmetry < 2);
        pruned = prepare_branch(i, current, it, nb, cdist, cTheta, sTheta, U, r1, r3, X2, space.v, S2, space.op, space.info);

        if (!pruned)
        {
          // Assign X2
          numChildren++;
        }
      }
      first = true;
    }
  }

  // Return the next node to look into
  SearchNode next() override
  {
    SearchNode nextNode;
    nextNode.i = i + 1;
    nextNode.n_vertices = n_vertices;
    if (!firstIsPruned && first)
    {
      nextNode.X = X1;
      nextNode.S = S1;
      first = false;
    }
    else
    {
      nextNode.X = X2;
      nextNode.S = S2;
    }
    return nextNode;
  }
};

int hpx_main(hpx::program_options::variables_map &opts)
{
  int i, n;
  VERTEX *v;
  double **X;
  SEARCH S;
  OPTION op;
  INFORMATION info;

  std::cout << "YewPar DDGP solver based on MDjeep" << std::endl;

  auto inputFile = opts["mdfile"].as<std::string>();

  Config config = mdjeep_main(inputFile);
  if (config.error)
    return hpx::finalize();
  n = config.n;
  m = config.m;
  v = config.v;
  X = config.X;
  S = config.S;
  op = config.op;
  info = *(config.info);

  i = prepare_bp(v, X, S, op, &info);
  // Ignoring bp_exact for now

  /*
    Main body
  */

  SearchNode root;
  root.i = i;
  root.n_vertices = n;
  root.X = X;
  root.S = S;
  SearchSpace searchS = {false, m, v, op, &info};

  fprintf(stderr, "mdjeep: bp is exploring the search tree ... ");
  if (op.monitor)
  {
    fprintf(stderr, "layer ");
    for (i = 0; i < info.ndigits; i++)
      fprintf(stderr, " ");
  };
  auto start_time = std::chrono::steady_clock::now();

  auto skeletonType = opts["skeleton"].as<std::string>();

  std::uint64_t count;
  if (skeletonType == "seq")
  {
    YewPar::Skeletons::API::Params<> searchParameters;
    count = YewPar::Skeletons::Seq<GenNode,
                                   YewPar::Skeletons::API::Enumeration,
                                   YewPar::Skeletons::API::Enumerator<CountSols>,
                                   YewPar::Skeletons::API::DepthLimited>::search(searchS, root, searchParameters);
  }
  else
  {
    hpx::cout << "Invalid skeleton type option. Only seq implemented so far." << std::endl;
    hpx::finalize();
    return EXIT_FAILURE;
  }

  auto overall_time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start_time);

  fprintf(stderr, "\n");
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