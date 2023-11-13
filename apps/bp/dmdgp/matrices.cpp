#include "bp.hpp"

// Function directly adapted from MD-jeep
std::array<double, 12> matrixProd(const std::array<double, 12> &qim1, const std::array<double, 12> &qi)
{
    std::array<double, 12> Qout;

    Qout[0] = qim1[0] * qi[0] + qim1[1] * qi[4] + qim1[2] * qi[8];
    Qout[1] = qim1[0] * qi[1] + qim1[1] * qi[5] + qim1[2] * qi[9];
    Qout[2] = qim1[0] * qi[2] + qim1[1] * qi[6] + qim1[2] * qi[10];
    Qout[3] = qim1[0] * qi[3] + qim1[1] * qi[7] + qim1[2] * qi[11] + qim1[3];

    Qout[4] = qim1[4] * qi[0] + qim1[5] * qi[4] + qim1[6] * qi[8];
    Qout[5] = qim1[4] * qi[1] + qim1[5] * qi[5] + qim1[6] * qi[9];
    Qout[6] = qim1[4] * qi[2] + qim1[5] * qi[6] + qim1[6] * qi[10];
    Qout[7] = qim1[4] * qi[3] + qim1[5] * qi[7] + qim1[6] * qi[11] + qim1[7];

    Qout[8] = qim1[8] * qi[0] + qim1[9] * qi[4] + qim1[10] * qi[8];
    Qout[9] = qim1[8] * qi[1] + qim1[9] * qi[5] + qim1[10] * qi[9];
    Qout[10] = qim1[8] * qi[2] + qim1[9] * qi[6] + qim1[10] * qi[10];
    Qout[11] = qim1[8] * qi[3] + qim1[9] * qi[7] + qim1[10] * qi[11] + qim1[11];

    return Qout;
}