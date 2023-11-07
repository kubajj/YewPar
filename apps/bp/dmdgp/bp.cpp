#include "bp.hpp"

DMDGPSol placeFirstThreeVertices(const std::map<std::pair<int, int>, DataRecord *> &references, const std::map<std::pair<int, int>, double> &thetaMap)
{
    DMDGPSol sol;

    // Matrices B1..B3
    // Qi = B4 \cdots Bi
    // last row 0,0,0,1 is never used, so we don't store it
    double b1[12], b2[12], b3[12];

    // First Vertex (0, 0, 0)
    DMDGPVertexPosition vertex1;
    vertex1.x = 0.0;
    vertex1.y = 0.0;
    vertex1.z = 0.0;
    // vertex1.atom = find the name of the atom;
    // vertex1.amino = find the name of the amino acid;
    sol.vertices.push_back(vertex1);
    b1 = [ 1, 0, 0, 0,
           0, 1, 0, 0,
           0, 0, 1, 0 ];

    // Second Vertex
    DMDGPVertexPosition vertex2;
    auto iref12 = references.find({1, 2});
    double d12;
    if (iref12 != references.end())
    {
        d12 = iref12->second->lb;
    }
    // else
    // {
    //     std::cerr << "Something went wrong when positioning second vertex" << std::endl;
    //     return NULL;
    // }
    vertex2.x = -d12;
    vertex2.y = 0.0;
    vertex2.z = 0.0;
    sol.vertices.push_back(vertex2);
    b2 = [ -1, 0, 0, -d12,
           0, 1, -1, 0,
           0, 0, 1, 0 ];

    // Third Vertex
    DMDGPVertexPosition vertex3;
    auto iref23 = references.find({2, 3});
    auto itheta13 = thetaMap.find({1, 3});
    double d23, theta13;
    // Cosine and sine of theta13
    double ct13, st13;
    if (it12 != references.end() && it23 != references.end() && it13 != thetaMap.end())
    {
        d23 = iref23->second->lb;
        theta13 = itheta13->second;
    }
    // else
    // {
    //     std::cerr << "Something went wrong when positioning third vertex" << std::endl;
    //     return NULL;
    // }
    ct13 = std::cos(theta13);
    st13 = std::sin(theta13);
    vertex3.x = -d12 + d23 * ct13;
    vertex3.y = d23 * st13;
    vertex3.z = 0.0;
    sol.vertices.push_back(vertex3);
    b3 = [ -ct13, -st13, 0, -d23 * ct13,
           st13, -ct13, 0, d23 * st13,
           0, 0, 1, 0 ];

    return sol;
}