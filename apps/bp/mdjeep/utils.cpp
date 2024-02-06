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

SEARCH copySearch(SEARCH S, int n, int m)
{
    SEARCH newS;
    newS.sym = (bool *)calloc(n, sizeof(bool));
    newS.refs = (triplet *)calloc(n, sizeof(triplet));
    S.pX = allocateMatrix(3, n);
    S.lX = allocateMatrix(3, n);
    S.uX = allocateMatrix(3, n);
    S.y = allocateVector(m);
    S.gy = allocateVector(m);
    S.sy = allocateVector(m);
    S.yp = allocateVector(m);
    S.gyp = allocateVector(m);
    S.gX = allocateMatrix(3, n);
    S.sX = allocateMatrix(3, n);
    S.Xp = allocateMatrix(3, n);
    S.gXp = allocateMatrix(3, n);
    S.DX = allocateMatrix(3, n);
    S.YX = allocateMatrix(3, n);
    S.ZX = allocateMatrix(3, n);
    S.Dy = allocateVector(m);
    S.Yy = allocateVector(m);
    S.Zy = allocateVector(m);
    S.memory = allocateVector(n);
    for (int i = 0; i < n; i++)
    {
        newS.sym[i] = S.sym[i];
        newS.refs[i] = S.refs[i];
        newS.memory[i] = S.memory[i];
        for (int j = 0; j < 3; j++)
        {
            newS.pX[j][i] = S.pX[j][i];
            newS.lX[j][i] = S.lX[j][i];
            newS.uX[j][i] = S.uX[j][i];
            newS.gX[j][i] = S.gX[j][i];
            newS.sX[j][i] = S.sX[j][i];
            newS.Xp[j][i] = S.Xp[j][i];
            newS.gXp[j][i] = S.gXp[j][i];
            newS.DX[j][i] = S.DX[j][i];
            newS.YX[j][i] = S.YX[j][i];
            newS.DX[j][i] = S.DX[j][i];
        }
    }
    for (int k = 0; k < m; k++)
    {
        newS.y[k] = S.y[k];
        newS.gy[k] = S.gy[k];
        newS.sy[k] = S.sy[k];
        newS.yp[k] = S.yp[k];
        newS.gyp[k] = S.gyp[k];
        newS.Dy[k] = S.Dy[k];
        newS.Yy[k] = S.Yy[k];
        newS.Zy[k] = S.Zy[k];
    }
    newS.pi = S.pi;
    return newS;
}