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

// Calculate the angle between atoms at indices i-3, i-2, i-1, and i using cosine rule
double calculateOmega(const std::map<std::pair<int, int>, double> &lbMap, int i)
{
    double a, b, c, e, f, cosOmega;
    int im3 = i - 3;
    int im2 = i - 2;
    int im1 = i - 1;

    // Check if the lb values exist in the map
    auto it1 = lbMap.find({im3, im2});
    auto it2 = lbMap.find({im3, im1});
    auto it3 = lbMap.find({im3, i});
    auto it4 = lbMap.find({im2, im1});
    auto it5 = lbMap.find({im2, i});
    auto it6 = lbMap.find({im1, i});

    if (it1 == lbMap.end() || it2 == lbMap.end() || it3 == lbMap.end() || it4 == lbMap.end() || it5 == lbMap.end() || it6 == lbMap.end())
    {
        throw std::runtime_error("Lower bounds not found for the quadruplet.");
    }

    // Retrieve the lower bounds (lb) values
    double d_im3_im2 = it1->second;
    double d_im3_im1 = it2->second;
    double d_im3_i = it3->second;
    double d_im2_im1 = it4->second;
    double d_im2_i = it5->second;
    double d_im1_i = it6->second;

    a = (d_im3_im2 * d_im3_im2 + d_im2_i * d_im2_i - d_im3_i * d_im3_i) / (2.0 * d_im3_im2 * d_im2_i);
    b = (d_im2_i * d_im2_i + d_im2_im1 * d_im2_im1 - d_im1_i * d_im1_i) / (2.0 * d_im2_i * d_im2_im1);
    c = (d_im3_im2 * d_im3_im2 + d_im2_im1 * d_im2_im1 - d_im3_im1 * d_im3_im1) / (2.0 * d_im3_im2 * d_im2_im1);
    e = 1.0 - b * b;
    f = 1.0 - c * c;
    if (e < 0.0 || f < 0.0)
    {
        throw std::runtime_error("Something went wrong during the computation of the cosine.");
    }
    e = sqrt(e);
    f = sqrt(f);
    cosOmega = (a - b * c) / (e * f);

    // Calculate and return the angle in radians
    return std::acos(cosOmega);
}

// Function to calculate and return a map of angles for vertices i from 3 to n
void calculateAnglesForVertices(
    const std::map<std::pair<int, int>, double> &lbMap, int n,
    std::map<std::pair<int, int>, double> &thetaMap,
    std::map<std::pair<int, int>, double> &omegaMap)
{

    try
    {
        double angle = calculateTheta(records, 3);
        // Store the angle in the angleMap with the key (im2, i)
        thetaMap[{1, 3}] = angle;
    }
    catch (const std::runtime_error &e)
    {
        std::cout << e.what() << std::endl;
    }

    for (int i = 4; i <= n; i++)
    {
        int im3 = i - 3;
        int im2 = i - 2;

        try
        {
            double theta = calculateTheta(records, i);
            // Store the angle in the angleMap with the key (im2, i)
            thetaMap[{im2, i}] = theta;
            double omega = calculateOmega(records, i);
            omegaMap[{im3, i}] = omega;
        }
        catch (const std::runtime_error &e)
        {
            std::cout << e.what() << std::endl;
        }
    }
}