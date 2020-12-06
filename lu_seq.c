
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

double **make2dmatrix(long n)
{
	long i;
	double **m;
	m = (double **)malloc(n * sizeof(double *));
	for (i = 0; i < n; i++)
		m[i] = (double *)malloc(n * sizeof(double));
	return m;
}

void printmatrix(double **A, long n)
{
	long i, j;
	for (i = 0; i < n; i++)
	{
		for (j = 0; j < n; j++)
			printf("%.2lf \t\t", A[i][j]);
		printf("\n");
	}
}

void free2dmatrix(double **M, long n)
{
	long i;
	if (!M)
		return;
	for (i = 0; i < n; i++)
		free(M[i]);
	free(M);
}

void lu_fat(double **A, long n)
{
	double **up;
	up = (double **)calloc(sizeof(double *), n);
	for (int i = 0; i < n; i++)
	{
		up[i] = (double *)calloc(sizeof(double), n);
		if (up[i] == NULL)
		{
			printf("CANNOT ALLOCATE, OUT OF MEMORY \n");
			return;
		}
	}

	double **lo;
	lo = (double **)calloc(sizeof(double *), n);
	for (int i = 0; i < n; i++)
	{
		lo[i] = (double *)calloc(sizeof(double), n);
		if (lo[i] == NULL)
		{
			printf("CANNOT ALLOCATE, OUT OF MEMORY \n");
			return;
		}
	}

	for (int i = 0; i < n; i++)
	{
		for (int k = i; k < n; k++)
		{
			double sum = 0;
			for (int j = 0; j < i; j++)
			{
				sum += (lo[i][j] * up[j][k]);
			}
			up[i][k] = A[i][k] - sum;
		}


		for (int k = i; k < n; k++)
		{
			if (i == k)
			{
				lo[i][i] = 1;
			}
			else
			{
				double sum = 0;
				for (int j = 0; j < i; j++)
				{
					sum += (lo[k][j] * up[j][i]);
				}
				lo[k][i] = (A[k][i] - sum) / up[i][i];
			}
		}
	}

	printf("LOWER\n");
	printmatrix(lo, n);
	printf("_________________\n");
	printf("UPPER\n");
	printmatrix(up, n);
}

void initialize(double **A, long n)
{
	time_t t;
	srand((unsigned)time(&t));
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i][j] = (double) (rand() % 100);
		}
	}
}
double **getMatrix(long size)
{
	double **m = make2dmatrix(size);
	initialize(m, size);
	return m;
}

int main()
{
	long matrix_size;
	printf("Enter the size of matrix (N x N) where N = ");
	scanf("%lu", &matrix_size);
	clock_t begin, end;
	double time_spent;

	double **matrix = getMatrix(matrix_size);
	printmatrix(matrix, matrix_size);
	begin = clock();
	lu_fat(matrix, matrix_size);
	end = clock();
	time_spent = ((double)(end - begin)) / CLOCKS_PER_SEC;
	free2dmatrix(matrix, matrix_size);
}