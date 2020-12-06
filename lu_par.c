
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
#include <upc_relaxed.h>
#include <upc.h>

shared[] double *m;

void lu_fat(shared[] double *A, long n)
{
	if (MYTHREAD == 0)
	{
		for (int i = 0; i < n; i++)
		{
			for (int j = i * n; j < i * n + n; j++)
			{
				printf("%.2f \t", A[j]);
			}
			printf("\n");
		}
		printf("-----------------------------------------\n");
	}

	int rows = n / THREADS;
	int min = MYTHREAD * rows;
	int max = min + rows - 1;

	if (MYTHREAD == (THREADS - 1) && (n - (max + 1)) > 0)
	{
		max = n - 1;
	}

	upc_forall(int k = 0; k < n; k++; continue)
	{
		if (k >= min && k <= max)
		{
			upc_forall(int j = k + 1; j < n; j++; continue)
			{
				A[(n * j) + k] = A[(n * j) + k] / A[(k * n) + k];
			}
		}
		upc_barrier;

		upc_forall(int i = (((k + 1) > min) ? (k + 1) : min); i <= max; i++; continue)
		{
			upc_forall(int j = k + 1; j < n; j++; continue)
			{

				A[(i * n) + j] = A[(i * n) + j] - (A[(i * n) + k] * A[(k * n) + j]);
			}
		}
		upc_barrier;
	}
}

void printmatrix(double **A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			printf("%.2f \t", A[i][j]);
		}
		printf("\n");
	}
}

void getMatrix(long size)
{
	m = upc_all_alloc(THREADS, size * size * sizeof(double));
	time_t t;
	srand((unsigned)time(&t));
	for (int i = 0; i < size * size; i++)
	{
		m[i] = rand() % 50;
	}
}

void parallelMultiply(double **matrixA, double **matrixB, double **matrixC, long n)
{
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			for (int k = 0; k < n; k++)
			{
				matrixC[i][j] += matrixA[i][k] * matrixB[k][j];
			}
		}
	}
}

double **make2dmatrix(long n)
{
	long i;
	double **m;
	m = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		m[i] = (double *)malloc(n * sizeof(double));
	return m;
}

void factorize(shared[] double *A, long n)
{
	long i, j;
	double **u = make2dmatrix(n);
	double **l = make2dmatrix(n);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			u[i][j] = A[(i * n) + j];
			if (i == j)
				l[i][j] = 1;
			else
				l[i][j] = 0;
		}
	}
	for (int i = 1; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			l[i][j] = A[(i * n) + j];
			u[i][j] = 0;
		}
	}

	printf("\nL Matrix :\n\n");
	printmatrix(l, n);
	printf("\nU Matrix :\n\n");
	printmatrix(u, n);
	double **a1 = make2dmatrix(n);
	parallelMultiply(l, u, a1, n);
	printf("\nA Matrix after Matrix Mutliplication :\n\n");
	printmatrix(a1, n);
}

int main()
{
	long matrix_size;
	matrix_size = 10;
	clock_t begin, end;
	double time_spent;
	getMatrix(matrix_size);
	begin = clock();
	lu_fat(m, matrix_size);
	end = clock();
	time_spent = ((double)(end - begin)) / CLOCKS_PER_SEC;
	if (MYTHREAD == 0)
	{
		factorize(m, matrix_size);
	}
}
