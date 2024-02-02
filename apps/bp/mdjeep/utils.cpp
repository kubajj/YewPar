#include "bp.hpp"

double **copy(double **matrix, int rows, int cols)
{
    double **copy = new double *[rows];
    for (int i = 0; i < rows; ++i)
    {
        copy[i] = new double[cols];
    }

    // Copy the elements from the original matrix to the new one
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < cols; ++j)
        {
            copy[i][j] = matrix[i][j];
        }
    }

    return copy;
}