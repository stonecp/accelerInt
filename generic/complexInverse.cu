#include "header.cuh"
#include "solver_props.cuh"
#include <cuComplex.h>

/* An implementation of complexLU factorization and inversion
   based on BLAS/LAPACK routines for CUDA

   Author: Nick Curtis
*/

__device__
int getComplexMax (const int n, const cuDoubleComplex * __restrict__ complexArr) {
    
    int maxInd = 0;
    if (n == 1)
        return maxInd;
    
    double maxVal = cuCabs(complexArr[INDEX(0)]);
    for (int i = 1; i < n; ++i) {
        if (cuCabs(complexArr[INDEX(i)]) > maxVal) {
            maxInd = i;
            maxVal = cuCabs(complexArr[INDEX(i)]);
        }
    }
    
    return maxInd;
}

///////////////////////////////////////////////////////////

__device__
void scaleComplex (const int n, const cuDoubleComplex val, cuDoubleComplex* arrX) {
    
    for (int i = 0; i < n; ++i) {
        arrX[INDEX(i)] = cuCmul(arrX[INDEX(i)], val);
    }
    
}

///////////////////////////////////////////////////////////

__device__
void swapComplex (const int n, cuDoubleComplex* arrX, const int incX, cuDoubleComplex* arrY, const int incY) {
    
    int ix = 0;
    int iy = 0;
    
    for (int i = 0; i < n; ++i) {
        cuDoubleComplex temp = arrX[INDEX(ix)];
        arrX[INDEX(ix)] = arrY[INDEX(iy)];
        arrY[INDEX(iy)] = temp;
        ix += incX;
        iy += incY;
    }
    
}

///////////////////////////////////////////////////////////

__device__
void complexGERU (const int n, const cuDoubleComplex alpha, const cuDoubleComplex* arrX,
                                    const cuDoubleComplex* arrY, const int incY, cuDoubleComplex* A, const int lda) {
    
    for (int j = 0; j < n; ++j) {
        if (cuCabs(arrY[INDEX(j * incY)]) > 0.0) {
      
            cuDoubleComplex temp = cuCmul(alpha, arrY[INDEX(j * incY)]);
      
            for (int i = 0; i < n; ++i) {
                A[INDEX(i + (lda * j))] = cuCfma(arrX[INDEX(i)], temp, A[INDEX(i + (lda * j))]);
            }
      
        }    
    }
    
}

///////////////////////////////////////////////////////////

__device__
void multiplyComplexUpperMV (const int n, cuDoubleComplex* x, const int lda, const cuDoubleComplex* A) {
    
    for (int j = 0; j < n; ++j) {
        if (cuCabs(x[j]) > 0.0) {
            cuDoubleComplex temp = x[INDEX(j)];
            for (int i = 0; i < j; ++i) {
                //x[i] += temp * A[i + (lda * j)];
                x[INDEX(i)] = cuCfma(temp, A[INDEX(i + (lda * j))], x[INDEX(i)]);
            }
            //x[j] *= A[j + (lda * j)];
            x[INDEX(j)] = cuCmul(x[INDEX(j)], A[INDEX(j + (lda * j))]);
        }
    }
    
}

///////////////////////////////////////////////////////////

__device__
void complexGEMV (const int m, const int n, const int lda, const cuDoubleComplex alpha, const cuDoubleComplex* A,
                                    const cuDoubleComplex* arrX, cuDoubleComplex* arrY) {
    
    // first: y = beta*y
    // beta = 1, so nothing
    
    // second: y = alpha*A*x + y
    
    for (int j = 0; j < n - 1; ++j) {

        if (cuCabs(arrX[j]) > 0.0) {
            cuDoubleComplex temp = cuCmul(alpha, arrX[INDEX(j)]);
            for (int i = 0; i < m; ++i) {
                //arrY[i] += temp * A[i + (m * j)];
                arrY[INDEX(i)] = cuCfma(temp, A[INDEX(i + (lda * j))], arrY[INDEX(i)]);
            }
        }
    }
    
}

///////////////////////////////////////////////////////////

__device__
void getComplexLU (const int n, const int LDA, cuDoubleComplex* __restrict__ A,
                    int* __restrict__ indPivot, int* __restrict__ info) {
    
    for (int j = 0; j < n; ++j) {
        
        // find pivot and test for singularity
        
        int jp = j + getComplexMax (n - j, &A[GRID_DIM * (j + (LDA * j))]);
        indPivot[INDEX(j)] = jp;

        if (cuCabs(A[INDEX(jp + (LDA * j))]) > 0.0) {
            
            // apply interchange to columns 1:n-1
            if (jp != j)
                swapComplex (n, &A[GRID_DIM * (j)], LDA, &A[GRID_DIM * (jp)], LDA);
            
            // compute elements j+1:m-1 of the jth column
            
            if (j < NSP - 1)
                scaleComplex (n - j - 1, cuCdiv(make_cuDoubleComplex(1.0, 0.0), A[INDEX(j + (LDA * j))]), &A[GRID_DIM * (j + 1 + (LDA * j))]);
            
        } else if (*info == 0) {
            *info = j;
            break;
        }
        
        // update trailing submatrix
        if (j < n - 1)
            complexGERU (n - j - 1, make_cuDoubleComplex(-1.0, 0.0), &A[GRID_DIM * (j + 1 + (LDA * j))], &A[GRID_DIM * (j + LDA * (j + 1))], LDA, &A[GRID_DIM * (j + 1 + LDA * (j + 1))], LDA);
        
    }
}

__device__
void getComplexInverseLU (const int n, const int LDA, cuDoubleComplex* __restrict__ A,
                            const int* __restrict__ indPivot,
                            cuDoubleComplex* __restrict__ work) {
    
    // form inv(U)
    for (int j = 0; j < n; ++j) {
        A[j + (LDA * j)] = cuCdiv(make_cuDoubleComplex(1.0, 0.0), A[j + (LDA * j)]);
        cuDoubleComplex Ajj = cuCmul(make_cuDoubleComplex(-1.0, 0.0), A[j + (LDA * j)]);
        
        // compute elements 0:j-1 of jth column
        multiplyComplexUpperMV (j, &A[LDA * j], LDA, A);
        
        // scale
        scaleComplex (j, Ajj, &A[LDA * j]);
    }
    
    // solve equation inv(A)*L = inv(U) for inv(A)
    
    for (int j = n - 1; j >= 0; --j) {
        
        // copy current column of L to work and replace with 0.0s
        for (int i = j + 1; i < n; ++i) {
            work[INDEX(i)] = A[INDEX(i + (LDA * j))];
            A[INDEX(i + (LDA * j))] = make_cuDoubleComplex(0.0, 0.0);
        }
        
        // compute current column of inv(A)
        if (j < n - 1)
            complexGEMV (n, n - j, LDA, make_cuDoubleComplex(-1.0, 0.0), &A[GRID_DIM * (LDA * (j + 1))], &work[GRID_DIM * (j + 1)], &A[GRID_DIM * (LDA * j)]);
        
    }
    
    // apply column interchanges
    
    for (int j = n - 2; j >= 0; --j) {
    
        if (indPivot[INDEX(j)] != j)
            swapComplex (n, &A[GRID_DIM * (LDA * j)], 1, &A[GRID_DIM * (LDA * indPivot[INDEX(j)])], 1);
    }
}

__device__
void getComplexInverse (const int n, const int LDA, cuDoubleComplex* __restrict__ A,
                            int* __restrict__ ipiv, int* __restrict__ info,
                            cuDoubleComplex* __restrict__ work) {

    // first get LU factorization
    getComplexLU (n, LDA, A, ipiv, info);

    // check for successful exit
    if (*info != 0) {
        return;
    }

    // now get inverse
    getComplexInverseLU (n, LDA, A, ipiv, work);
}

//Matrix Algorithms: Volume 1: Basic Decompositions
//By G. W. Stewart
__device__
void getHessenbergLU(const int n, const int LDA, cuDoubleComplex* A, int* __restrict__ indPivot, int* __restrict__ info)
{
    int last_free = 0;
    for (int i = 0; i < n - 1; i ++)
    {
        if (cuCabs(A[INDEX(i * LDA + i)]) < cuCabs(A[i * LDA + i + 1]))
        {
            //swap rows
            swapComplex(n - last_free, &A[GRID_DIM * (last_free * LDA + i)], LDA, &A[GRID_DIM * (last_free * LDA + i + 1)], LDA);
            indPivot[INDEX(i)] = i + 1;
        }
        else
        {
            indPivot[INDEX(i)] = i;
            last_free = i;
        }
        if (cuCabs(A[INDEX(i * LDA + i)]) > 0.0)
        {
            cuDoubleComplex tau = cuCdiv(A[INDEX(i * LDA + i + 1)], A[INDEX(i * LDA + i)]);
            for (int j = i + 1; j < n; j++)
            {
                A[INDEX(j * LDA + i + 1)] = cuCsub(A[INDEX(j * LDA + i + 1)], cuCmul(tau, A[INDEX(j * LDA + i)]));
            }
            A[INDEX(i * LDA + i + 1)] = tau;
        }
        else
        {
            *info = i;
            return;
        }
    }
    //last index is not pivoted
    indPivot[INDEX(n - 1)] = n - 1;
}

__device__
void getComplexInverseHessenberg (const int n, const int LDA, cuDoubleComplex* __restrict__ A,
                                    int* __restrict__ ipiv, int* __restrict__ info,
                                    cuDoubleComplex* __restrict__ work)
{
    // first get LU factorization
    getHessenbergLU (n, LDA, A, ipiv, info);

    if (*info != 0)
        return;

    // now get inverse
    getComplexInverseLU (n, LDA, A, ipiv, work);
}