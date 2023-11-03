#include "bp.hpp"

ParsedData parseFile(const std::string &filename)
{
    ParsedData data;
    std::ifstream file(filename);
    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file");
    }

    std::map<std::string, std::string> expectedKeys;

    expectedKeys["instance"] = "";
    expectedKeys["with file"] = "";
    // expectedKeys["with format"] = "";
    // expectedKeys["with separator"] = "";
    expectedKeys["method"] = "";
    expectedKeys["with resolution"] = "";
    expectedKeys["with tolerance"] = "";
    expectedKeys["with maxtime"] = "";

    std::string line;
    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            // Skip empty lines and comments
            continue;
        }

        size_t colonPos = line.find(':');
        if (colonPos != std::string::npos)
        {
            std::string key = line.substr(0, colonPos);
            std::string value = line.substr(colonPos + 1);
            value.erase(0, value.find_first_not_of(" \t")); // Trim leading spaces and tabs

            if (expectedKeys.find(key) != expectedKeys.end())
            {
                expectedKeys[key] = value;
            }
        }
    }

    data.file = expectedKeys["with file"];
    // data.format = expectedKeys["with format"];
    // data.separator = expectedKeys["with separator"];
    data.methodParams.resolution = std::stod(expectedKeys["with resolution"]);
    data.methodParams.tolerance = std::stod(expectedKeys["with tolerance"]);
    data.methodParams.maxtime = std::stoi(expectedKeys["with maxtime"]);

    return data;
}

std::vector<DataRecord> readDataFile(const std::string &filename, int &maxId)
{
    std::vector<DataRecord> records;
    std::ifstream file(filename);

    if (!file.is_open())
    {
        throw std::runtime_error("Failed to open file");
    }

    std::string line;
    maxId = 0; // Initialize maxId

    while (std::getline(file, line))
    {
        if (line.empty() || line[0] == '#')
        {
            continue; // Skip empty lines and comments
        }

        DataRecord record;
        std::istringstream tokenStream(line);

        tokenStream >> record.Id1 >> record.Id2 >> record.groupId1 >> record.groupId2 >>
            record.lb >> record.ub >> record.Name1 >> record.Name2 >> record.groupName1 >> record.groupName2;

        records.push_back(record);

        // Update maxId with the maximum of Id1 and Id2
        maxId = std::max({maxId, record.Id1, record.Id2});
    }

    return records;
}

std::map<std::pair<int, int>, DataRecord *> createDataRecordMap(const std::vector<DataRecord> &records)
{
    std::map<std::pair<int, int>, DataRecord *> recordMap;

    for (auto &record : records)
    {
        // Create a key as a pair of Id1 and Id2
        std::pair<int, int> key(record.Id1, record.Id2);

        // Associate the key with a pointer to the data record
        recordMap[key] = &record;
    }

    return recordMap;
}