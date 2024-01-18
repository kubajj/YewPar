#include "bp.hpp"

bool pruningTest(int vertexId, const DMDGPMaps &maps, DMDGPSol &sol, DMDGPVertexPosition currentPosition)
{
    int j = vertexId;
    bool prune = false;
    double eps, dist, diff;
    // Epsilon for precision
    eps = 0.001;
    // Note(kubajj): Vertices are numbered 1..n,
    // but sol.vertices starts at 0
    for (int i = 1; i < j && !prune; i++)
    {
        auto idistance = maps.distanceMap.find({i, vertexId});
        // There is a distance between i and current vertex
        if (idistance != maps.distanceMap.end())
        {
            dist = calculateDistance(sol.vertices[i - 1], currentPosition);
            diff = dist - idistance->second;
            if (diff < 0.0)
                diff = -diff;
            if (diff > eps)
            {
                prune = true;
            }
            // Can calculate partial lde here
        }
    }
    return prune;
}