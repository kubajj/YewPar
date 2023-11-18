#include "bp.hpp"

DMDGPSol placeFirstThreeVertices(const std::map<std::pair<int, int>, double> &distanceMap, const std::map<std::pair<int, int>, double> &cosThetaMap, std::array<double, 12> &q3)
{
  DMDGPSol sol;

  // Matrices B1..B3
  // Qi = B4 \cdots Bi
  // last row 0,0,0,1 is never used, so we don't store it

  // First Vertex (0, 0, 0)
  DMDGPVertexPosition vertex1;
  vertex1.x = 0.0;
  vertex1.y = 0.0;
  vertex1.z = 0.0;
  // vertex1.atom = find the name of the atom;
  // vertex1.amino = find the name of the amino acid;
  sol.vertices.push_back(vertex1);
  std::array<double, 12> b1 = {1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0};

  // Second Vertex
  DMDGPVertexPosition vertex2;
  auto iref12 = distanceMap.find({1, 2});
  double d12;
  d12 = iref12->second;
  vertex2.x = -d12;
  vertex2.y = 0.0;
  vertex2.z = 0.0;
  sol.vertices.push_back(vertex2);
  std::array<double, 12> b2 = {-1, 0, 0, -d12, 0, 1, -1, 0, 0, 0, 1, 0};
  std::array<double, 12> q2 = matrixProd(b1, b2);

  // Third Vertex
  DMDGPVertexPosition vertex3;
  auto iref23 = distanceMap.find({2, 3});
  auto itheta13 = cosThetaMap.find({1, 3});
  double d23, theta13;
  // Cosine and sine of theta13
  double ct13, st13;
  d23 = iref23->second;
  ct13 = itheta13->second;
  st13 = sqrt(1.0 - ct13 * ct13);
  vertex3.x = -d12 + d23 * ct13;
  vertex3.y = d23 * st13;
  vertex3.z = 0.0;
  sol.vertices.push_back(vertex3);
  std::array<double, 12> b3 = {-ct13, -st13, 0, -d23 * ct13, st13, -ct13, 0, d23 * st13, 0, 0, 1, 0};
  // q3 = {-ct13, -st13, 0, -d23 * ct13,
  //       st13, -ct13, 0, d23 * st13,
  //       0, 0, 1, 0};
  q3 = matrixProd(q2, b3);

  return sol;
}

void calculateBis(int i, std::array<double, 12> bi1, std::array<double, 12> bi2, const DMDGPMaps maps)
{
  int im1 = i - 1;
  int im2 = i - 2;
  int im3 = i - 3;
  auto idistance = maps.distanceMap.find({im1, i});
  auto itheta = maps.cosThetaMap.find({im2, i});
  auto iomega = maps.cosOmegaMap.find({im3, i});
  double d, cosTheta, sinTheta, cosOmega, sinOmega;
  d = idistance->second;
  cosTheta = itheta->second;
  cosOmega = iomega->second;
  sinTheta = sqrt(1.0 - cosTheta * cosTheta);
  sinOmega = sqrt(1.0 - cosOmega * cosOmega);
  // clang-format off
  bi1 = {
    -cosTheta, -sinTheta, 0, -d*cosTheta,
    sinTheta*cosOmega, -cosTheta*cosOmega, -sinOmega, d*sinTheta*cosOmega,
    sinTheta*sinOmega, -cosTheta*sinOmega, cosOmega, d*sinTheta*sinOmega
  };
  bi2 = {
    -cosTheta, -sinTheta, 0, -d*cosTheta,
    sinTheta*cosOmega, -cosTheta*cosOmega, sinOmega, d*sinTheta*cosOmega,
    -sinTheta*sinOmega, cosTheta*sinOmega, -cosOmega, d*sinTheta*sinOmega
  };
  // clang-format on
}