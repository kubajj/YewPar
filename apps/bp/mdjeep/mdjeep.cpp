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
*****************************************************************************************************/

/*
    MDjeep main.c - Adapted for integration with YewPar
*/

#include "bp.hpp"

Config mdjeep_main(std::string inputFile)
{
    Config config;
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

    // verifying whether all distances are exact (and precise)
    if (m == mexact)
    {
        fprintf(stderr, "mdjeep: the instance contains only 'exact' distances\n");
        if (info.method == 0)
        {
            // we'll invoke bp_exact if at least 90% of the distances are "very precise"
            if (totalNumberOfPreciseDistances(n, v, 14) > 0.90 * mexact)
            {
                op.r = 0.0;
                info.exact = true;
                fprintf(stderr, "mdjeep: the resolution parameter and the refinement method have been disabled\n");
            };
        };
    };

    // error if no refinement method was specified for bp, and the instance does not contain enough precise distances
    if (!info.exact)
    {
        if (info.method == 0)
        {
            if (info.refinement == -1)
            {
                fprintf(stderr, "mdjeep: error: no refinement method specified for bp (instance contains interval distances)\n");
                free(v);
                return hpx::finalize();
            };
        };
    };

    // checking whether the first three instance vertices form a clique (prerequisite for bp)
    if (info.method == 0)
    {
        clique = initialClique(n, v, op.eps);
        if (!clique)
        {
            fprintf(stderr, "mdjeep: error: the first three vertices of the input instance do not form a clique\n");
            fprintf(stderr, "               the instance cannot be discretized\n");
            free(v);
            return hpx::finalize();
        };
    };

    // checking whether the input instance is discretizable (prerequisite for bp)
    if (info.method == 0)
    {
        i = isDDGP(n, v, op.eps, clique);
        if (i != 0)
        {
            fprintf(stderr, "mdjeep: error: the input instance is not discretizable\n");
            fprintf(stderr, "               not enough references for vertex %d (should have at least 3, at least 2 exact)\n", n0 + i);
            fprintf(stderr, "               stopping here... other necessary distances may be unavailable\n");
            free(v);
            return hpx::finalize();
        };
    };

    // if bp is selected, we know that the input instance is discretizable
    if (info.method == 0)
        fprintf(stderr, "mdjeep: the input instance is discretizable\n");

    // checking the consecutivity assumption (optional)
    if (info.method == 0)
    {
        if (info.exact || check_consec)
        {
            info.consec = isDMDGP(n, v, op.eps, true);
            fprintf(stderr, "mdjeep: the instance ");
            if (info.consec)
                fprintf(stderr, "satisfies ");
            else
                fprintf(stderr, "does not satisfy ");
            fprintf(stderr, "the consecutivity assumption\n");
        };
    };

    // preparing for calling bp method
    if (info.method == 0)
    {
        // memory allocation for S.refs
        S.refs = (triplet *)calloc(n, sizeof(triplet));

        // initializing all reference triplets to null
        for (i = 0; i < n; i++)
            S.refs[i] = nullTriplet();

        // the definition of the reference vertices depends on the presence of interval distances
        // for instances with exact distances only: the code below only verifies that the flattest triplet is not "too flat"
        smallsine = false;
        for (i = 3; i < n; i++)
        {
            if (info.exact || onlyPreciseDistances(v[i].ref, 14))
            {
                // triplet of references with exact distances
                cosine = 0.0;
                S.refs[i] = findReferencesExactCase(i, v, op.eps, &cosine);
                if (isNullTriplet(S.refs[i]))
                {
                    fprintf(stderr, "mdjeep: internal error: it was verified that the discretization assumptions were satisfied but they are actually not\n");
                    free(S.refs);
                    free(v);
                    return hpx::finalize();
                };
                if (cosine == 0.0)
                {
                    fprintf(stderr, "mdjeep: error: one triplet of reference vertices forms a flat angle; no alternative triplet available\n");
                    free(S.refs);
                    free(v);
                    return hpx::finalize();
                };
                if (fabs(sqrt(1.0 - cosine * cosine)) < op.eps)
                    smallsine = true;
            }
            else
            {
                // triplet of references with one interval distance
                S.refs[i] = findReferencesIntervalCase(i, v, op.eps);
                if (isNullTriplet(S.refs[i]))
                {
                    fprintf(stderr, "mdjeep: internal error: it was verified that the discretization assumptions were satisfied but they are actually not\n");
                    free(S.refs);
                    free(v);
                    return hpx::finalize();
                };
            };
        };
        if (smallsine)
        {
            fprintf(stderr, "mdjeep: WARNING: some triplets of reference vertices form a angle whose sine is very close to zero (tolerance is %g)\n", op.eps);
        };
    };

    // memory allocation for the solution X
    X = allocateMatrix(3, n);

    // preparing for calling spg method
    if (info.method == 1)
    {
        // reading the starting point from a text file (with predefined format)
        input = fopen(info.start, "r");
        if (input == NULL)
        {
            fprintf(stderr, "mdjeep: error while opening file containing starting point for spg\n");
            free(X);
            free(v);
            return hpx::finalize();
        };
        if (readStartingPoint(input, n, X) != n)
        {
            fprintf(stderr, "mdjeep: error while reading starting point for spg, it seems it doesnt contain the expected number of vertex positions\n");
            free(X);
            free(v);
            return hpx::finalize();
        };
    };

    // looking for symmetries (even when main method is spg)
    S.sym = (bool *)calloc(n, sizeof(bool));
    fprintf(stderr, "mdjeep: checking symmetries ... ");
    findSymmetries(n, v, S.sym);
    fprintf(stderr, "layers:");
    for (i = 0; i < n; i++)
        if (S.sym[i])
            fprintf(stderr, " %d", n0 + i);
    fprintf(stderr, "\n");

    // counting the maximum number of digits for monitor (optional)
    if (op.monitor)
    {
        if (info.method == 0)
            info.ndigits = numberOfDigits(n);
        else
            info.ndigits = numberOfDigits(op.maxit);
    };

    // setting up output filename (if necessary)
    if (op.print != 0)
        info.output = removExtension(info.filename);

    // message about additional memory allocation
    fprintf(stderr, "mdjeep: allocating memory ...");

    // memory allocation for the arrays in SEARCH (for both bp and spg)
    S.pX = allocateMatrix(3, n);
    S.lX = allocateMatrix(3, n);
    S.uX = allocateMatrix(3, n);
    S.y = allocateVector(m);
    S.gy = allocateVector(m);
    S.sy = allocateVector(m);
    S.yp = allocateVector(m);
    S.gyp = allocateVector(m);
    S.gX = allocateMatrix(3, n);
    S.sX = allocateMatrix(3, n);
    S.Xp = allocateMatrix(3, n);
    S.gXp = allocateMatrix(3, n);
    S.DX = allocateMatrix(3, n);
    S.YX = allocateMatrix(3, n);
    S.ZX = allocateMatrix(3, n);
    S.Dy = allocateVector(m);
    S.Yy = allocateVector(m);
    S.Zy = allocateVector(m);
    S.memory = allocateVector(n);
    fprintf(stderr, " done\n");

    // setting up value for pi
    S.pi = 3.14159265358979323846;

    // preparing for calling spg method (continues)
    if (info.method == 1)
    {
        // defining boxes necessary for spg
        for (i = 0; i < n; i++)
        {
            createBox(i, X, op.be, S.lX, S.uX);
            expandBounds(i, v, S.lX, S.uX, op.be, op.eps);
        };
    };
    config = config.createConfig(n, v, X, S, op, &info);
    return config;
}