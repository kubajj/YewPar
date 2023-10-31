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

ParsedData parseFile(const std::string &filename);
std::vector<DataRecord> readDataFile(const std::string &filename);