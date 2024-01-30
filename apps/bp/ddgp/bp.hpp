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

struct DDGPVertexPosition
{
    // char atom[3];
    // char amino[4];
    double x;
    double y;
    double z;
};

struct DDGPSol
{
    std::vector<DDGPVertexPosition> vertices;
};

struct DistanceMap
{
    std::map<std::pair<int, int>, std::pair<double, double>> boundsMap;

    double lowerBound(int key1, int key2) const
    {
        std::pair<int, int> keys(key1, key2);
        auto it = boundsMap.find(keys);

        if (it != boundsMap.end())
        {
            return it->second.first;
        }
        else
        {
            std::cerr << "Keys not found in the map." << std::endl;
            return NULL;
        }
    }

    double upperBound(int key1, int key2) const
    {
        std::pair<int, int> keys(key1, key2);
        auto it = boundsMap.find(keys);

        if (it != boundsMap.end())
        {
            return it->second.second;
        }
        else
        {
            std::cerr << "Keys not found in the map." << std::endl;
            return NULL;
        }
    }
};

struct VertexReferences
{
    std::map<int, std::vector<int>> references;

    std::vector<int> &getRefs(int key)
    {
        return references[key];
    }

    void insert(int key, int value)
    {
        references[key].push_back(value);
    }
};

struct DDGPMaps
{
    DistanceMap distanceMap;
    std::map<std::pair<int, int>, double> cosThetaMap;
    std::map<std::pair<int, int>, double> cosOmegaMap;
};

// readfile.cpp
ParsedData parseFile(const std::string &filename);
std::vector<DataRecord> readDataFile(const std::string &filename, int &maxId);
DistanceMap createDataRecordMap(std::vector<DataRecord> &records, VertexReferences &refs);