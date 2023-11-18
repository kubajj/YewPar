#ifndef _BP_
#define _BP_
#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
#include <array>
#include <map>
#include <chrono>
#include <memory>
#include <typeinfo>
#include <string>
#include <fstream>
#include <sstream>
#include <cmath>
#endif

struct MethodParameters
{
   double resolution;
   double tolerance;
   int maxtime;
};

struct ParsedData
{
   std::string file;
   // std::string format;
   // std::string separator;
   MethodParameters methodParams;
};

// Line in an instance file
struct DataRecord
{
   int Id1;
   int Id2;
   int groupId1;
   int groupId2;
   double lb;
   double ub;
   std::string Name1;
   std::string Name2;
   std::string groupName1;
   std::string groupName2;
};

struct DMDGPVertexPosition
{
   // char atom[3];
   // char amino[4];
   double x;
   double y;
   double z;
};

struct DMDGPSol
{
   std::vector<DMDGPVertexPosition> vertices;
};

struct DMDGPMaps
{
   std::map<std::pair<int, int>, double> distanceMap;
   std::map<std::pair<int, int>, double> cosThetaMap;
   std::map<std::pair<int, int>, double> cosOmegaMap;
   int number_of_vertices;
};

// readfile.cpp
ParsedData parseFile(const std::string &filename);
std::vector<DataRecord> readDataFile(const std::string &filename, int &maxId);
std::map<std::pair<int, int>, double> createDataRecordMap(std::vector<DataRecord> &records);

// vertex.cpp
void calculateAnglesForVertices(
    const std::map<std::pair<int, int>, double> &distanceMap, int n,
    std::map<std::pair<int, int>, double> &cosThetaMap,
    std::map<std::pair<int, int>, double> &cosOmegaMap);
double calculateDistance(DMDGPVertexPosition v1, DMDGPVertexPosition v2);

// bp.cpp
DMDGPSol placeFirstThreeVertices(
    const std::map<std::pair<int, int>, double> &distanceMap,
    const std::map<std::pair<int, int>, double> &cosThetaMap,
    std::array<double, 12> &q3);
void calculateBis(int i, std::array<double, 12> &bi1, std::array<double, 12> &bi2, const DMDGPMaps maps);

// pruningtest.cpp
bool pruningTest(int vertexId, const DMDGPMaps &maps, DMDGPSol &sol, DMDGPVertexPosition currentPosition);

// matrices.cpp
std::array<double, 12> matrixProd(const std::array<double, 12> &qim1, const std::array<double, 12> &qi);