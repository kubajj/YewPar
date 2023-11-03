#include "bp.hpp"

// Calculate the angle between consecutive triplets of atoms using cosine rule
double calculateTheta(const std::map<std::pair<int, int>, double> &records, int i)
{
    int im1 = i - 1;
    int im2 = i - 2;

    // Check if the lb values exist in the map
    auto it1 = records.find({im2, im1});
    auto it2 = records.find({im1, i});
    auto it3 = records.find({im2, i});

    if (it1 == records.end() || it2 == records.end() || it3 == records.end())
    {
        throw std::runtime_error("Lower bounds not found for the triplet.");
    }

    // Retrieve the lower bounds (lb) values
    // Note(kubajj): it->first is the key and it->second is the value
    double d_im2_im1 = it1->second->lb;
    double d_im1_i = it2->second->lb;
    double d_im2_i = it3->second->lb;

    // Compute the cosine of the angle
    // Using cosine rule:
    // cos angle = (a^2 + b^2 - c^2)/(2 * a * b)
    // The angle is between a and b
    double cosAngle = (d_im2_im1 * d_im2_im1 + d_im1_i * d_im1_i - d_im2_i * d_im2_i) / (2.0 * d_im2_im1 * d_im1_i);

    // Calculate and return the angle in radians
    return std::acos(cosAngle);
}

// Function to calculate and return a map of angles for vertices i from 3 to n
std::map<std::pair<int, int>, double> calculateThetasForVertices(const std::map<std::pair<int, int>, double> &records, int n)
{
    std::map<std::pair<int, int>, double> angleMap;

    for (int i = 3; i <= n; i++)
    {
        int im2 = i - 2;

        try
        {
            double angle = calculateTheta(records, i);
            // Store the angle in the angleMap with the key (im2, i)
            angleMap[{im2, i}] = angle;
        }
        catch (const std::runtime_error &e)
        {
            std::cout << e.what() << std::endl;
        }
    }

    return angleMap;
}