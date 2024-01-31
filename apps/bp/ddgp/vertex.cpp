#include "bp.hpp"

double calculateCosTheta(DistanceMap &distanceMap, int i, int im1, int im2)
{
    // Retrieve the lower bounds (lb) values
    // Note(kubajj): it->first is the key and it->second is the value
    double d_im2_im1 = distanceMap.lowerBound(im2, im1);
    double d_im1_i = distanceMap.lowerBound(im1, i);
    double d_im2_i = distanceMap.lowerBound(im2, i);

    // Compute the cosine of the angle
    // Using cosine rule:
    // cos angle = (a^2 + b^2 - c^2)/(2 * a * b)
    // The angle is between a and b
    double cosAngle = (d_im2_im1 * d_im2_im1 + d_im1_i * d_im1_i - d_im2_i * d_im2_i) / (2.0 * d_im2_im1 * d_im1_i);

    // Calculate and return the angle in radians
    return cosAngle;
}

double calculateCosOmega(DistanceMap &distanceMap, int i, int im1, int im2, int im3)
{
    double a, b, c, e, f, cosOmega;
    double d_im3_im2 = distanceMap.lowerBound(im3, im2);
    double d_im3_im1 = distanceMap.lowerBound(im3, im1);
    double d_im3_i = distanceMap.lowerBound(im3, i);
    double d_im2_im1 = distanceMap.lowerBound(im2, im1);
    double d_im2_i = distanceMap.lowerBound(im2, i);
    double d_im1_i = distanceMap.lowerBound(im1, i);
    a = (d_im3_im2 * d_im3_im2 + d_im2_i * d_im2_i - d_im3_i * d_im3_i) / (2.0 * d_im3_im2 * d_im2_i);
    b = (d_im2_i * d_im2_i + d_im2_im1 * d_im2_im1 - d_im1_i * d_im1_i) / (2.0 * d_im2_i * d_im2_im1);
    c = (d_im3_im2 * d_im3_im2 + d_im2_im1 * d_im2_im1 - d_im3_im1 * d_im3_im1) / (2.0 * d_im3_im2 * d_im2_im1);
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