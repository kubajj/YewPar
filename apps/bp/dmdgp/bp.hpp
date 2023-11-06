#ifndef _BP_
#define _BP_
#include <iostream>
#include <numeric>
#include <algorithm>
#include <vector>
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

struct DMDGPNode
{
   int id;
   DMDGPSol sol;

   int getObj() const
   {
      return id;
   }
};

// readfile.cpp
ParsedData parseFile(const std::string &filename);
std::vector<DataRecord> readDataFile(const std::string &filename, int &maxId);
std::map<std::pair<int, int>, DataRecord *> createDataRecordMap(std::vector<DataRecord> &records);

// vertex.cpp
void calculateAnglesForVertices(
    const std::map<std::pair<int, int>, DataRecord *> &references, int n,
    std::map<std::pair<int, int>, double> &thetaMap,
    std::map<std::pair<int, int>, double> &omegaMap);

// bp.cpp
DMDGPSol placeFirstThreeVertices(const std::map<std::pair<int, int>, DataRecord *> &references, const std::map<std::pair<int, int>, double> &thetaMap);

// pruningtest.cpp
bool pruningtest(int natoms, const std::vector<DataRecord> &records, std::map<std::pair<int, int>, DataRecord *> &references, DMDGPSol &sol);