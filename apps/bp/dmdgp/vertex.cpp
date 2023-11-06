#include "bp.hpp"

// Calculate the angle between consecutive triplets of atoms using cosine rule
double calculateCosTheta(const std::map<std::pair<int, int>, DataRecord *> &references, int i, int im1, int im2)
{

    // Check if the lb values exist in the map
    auto it1 = references.find({im2, im1});
    auto it2 = references.find({im1, i});
    auto it3 = references.find({im2, i});

    if (it1 == references.end() || it2 == references.end() || it3 == references.end())
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
    return cosAngle;
}

// Calculate the angle between atoms at indices i-3, i-2, i-1, and i using cosine rule
double calculateCosOmega(const std::map<std::pair<int, int>, DataRecord *> &references, int i, int im1, int im2, int im3)
{
    double a, b, c, e, f, cosOmega;

    a = calculateCosTheta(references, i, im2, im3);
    b = calculateCosTheta(references, im1, im2, i);
    c = calculateCosTheta(references, im1, im2, im3);
    //  a = (d_im3_im2 * d_im3_im2 + d_im2_i * d_im2_i - d_im3_i * d_im3_i) / (2.0 * d_im3_im2 * d_im2_i);
    //  b = (d_im2_i * d_im2_i + d_im2_im1 * d_im2_im1 - d_im1_i * d_im1_i) / (2.0 * d_im2_i * d_im2_im1);
    //  c = (d_im3_im2 * d_im3_im2 + d_im2_im1 * d_im2_im1 - d_im3_im1 * d_im3_im1) / (2.0 * d_im3_im2 * d_im2_im1);
    e = 1.0 - b * b;
    f = 1.0 - c * c;
    if (e < 0.0 || f < 0.0)
    {
        throw std::runtime_error("Something went wrong during the computation of the cosine of omega.");
    }
    e = sqrt(e);
    f = sqrt(f);
    cosOmega = (a - b * c) / (e * f);

    // Calculate and return the angle in radians
    return cosOmega;
}

// Function to calculate and return a map of angles for vertices i from 3 to n
void calculateAnglesForVertices(
    const std::map<std::pair<int, int>, DataRecord *> &references, int n,
    std::map<std::pair<int, int>, double> &thetaMap,
    std::map<std::pair<int, int>, double> &omegaMap)
{

    try
    {
        double angle = std::acos(calculateCosTheta(references, 3, 2, 1));
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
        int im1 = i - 1;

        try
        {
            double theta = std::acos(calculateCosTheta(references, i, im1, im2));
            // Store the angle in the angleMap with the key (im2, i)
            thetaMap[{im2, i}] = theta;
            double omega = std::acos(calculateCosOmega(references, i, im1, im2, im3));
            omegaMap[{im3, i}] = omega;
        }
        catch (const std::runtime_error &e)
        {
            std::cout << e.what() << std::endl;
        }
    }
}