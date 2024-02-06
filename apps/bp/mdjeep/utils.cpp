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
        if (S.sym[i] != NULL)
            newS.sym[i] = S.sym[i];
        else
            newS.sym[i] = NULL;
        if (!isNullTriplet(S.refs[i]))
            newS.refs[i] = S.refs[i];
        else
            newS.refs[i] = nullTriplet();
        if (S.memory[i] != NULL)
            newS.memory[i] = S.memory[i];
        else
            newS.memory[i] = NULL;
        for (int j = 0; j < 3; j++)
        {
            if (S.pX[j][i] != NULL)
                newS.pX[j][i] = S.pX[j][i];
            else
                newS.pX[j][i] = NULL;
            if (S.lX[j][i] != NULL)
                newS.lX[j][i] = S.lX[j][i];
            else
                newS.lX[j][i] = NULL;
            if (S.uX[j][i] != NULL)
                newS.uX[j][i] = S.uX[j][i];
            else
                newS.uX[j][i] = NULL;
            if (S.gX[j][i] != NULL)
                newS.gX[j][i] = S.gX[j][i];
            else
                newS.gX[j][i] = NULL;
            if (S.sX[j][i] != NULL)
                newS.sX[j][i] = S.sX[j][i];
            else
                newS.sX[j][i] = NULL;
            if (S.Xp[j][i] != NULL)
                newS.Xp[j][i] = S.Xp[j][i];
            else
                newS.Xp[j][i] = NULL;
            if (S.gXp[j][i] != NULL)
                newS.gXp[j][i] = S.gXp[j][i];
            else
                newS.gXp[j][i] = NULL;
            if (S.DX[j][i] != NULL)
                newS.DX[j][i] = S.DX[j][i];
            else
                newS.DX[j][i] = NULL;
            if (S.YX[j][i] != NULL)
                newS.YX[j][i] = S.YX[j][i];
            else
                newS.YX[j][i] = NULL;
            if (S.ZX[j][i] != NULL)
                newS.ZX[j][i] = S.ZX[j][i];
            else
                newS.ZX[j][i] = NULL;
        }
    }
    for (int k = 0; k < m; k++)
    {
        if (S.y[k] != NULL)
            newS.y[k] = S.y[k];
        else
            newS.y[k] = NULL;
        if (S.gy[k] != NULL)
            newS.gy[k] = S.gy[k];
        else
            newS.gy[k] = NULL;
        if (S.sy[k] != NULL)
            newS.sy[k] = S.sy[k];
        else
            newS.sy[k] = NULL;
        if (S.yp[k] != NULL)
            newS.yp[k] = S.yp[k];
        else
            newS.yp[k] = NULL;
        if (S.gyp[k] != NULL)
            newS.gyp[k] = S.gyp[k];
        else
            newS.gyp[k] = NULL;
        if (S.Dy[k] != NULL)
            newS.Dy[k] = S.Dy[k];
        else
            newS.Dy[k] = NULL;
        if (S.Yy[k] != NULL)
            newS.Yy[k] = S.Yy[k];
        else
            newS.Yy[k] = NULL;
        if (S.Zy[k] != NULL)
            newS.Zy[k] = S.Zy[k];
        else
            newS.Zy[k] = NULL;
    }
    newS.pi = S.pi;
    return newS;
}