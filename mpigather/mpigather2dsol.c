#include <stdio.h>
#include <mpi.h>

/*
 * This extends the 1D "cyclicgather" example to a 2D distributed array.
 * For example, consider a 4x4 matrix distribute across a 2x2 process grid:
 *
 *     Matrix             Process grid
 *
 *   0  1 |  2  3      
 *   4  5 |  6  7        rank 0 | rank 1
 *   -----+------        -------+-------
 *   8  9 | 10 11        rank 2 | rank 3
 *  12 13 | 14 15 
 *
 * Use a combination of MPI_Gatherv and a resized vector datatype to
 * reassemble this in the correct order on rank 0.
 */

#define M 12    // Global number of rows
#define N 10    // Global number of columns

#define ROWP 3  // Number of processes for row distribution
#define COLP 2  // Number of processes for column distribution

#define NPROC (ROWP*COLP)  // Total number of processes

#define MP (M/ROWP)  // Local number of rows
#define NP (N/COLP)  // Local number of columns

void matprint(int *mat, int m, int n);

void main(void)
{
  int i, j, iglobal, jglobal;
  int rank, size;
  int ipos, jpos;
  
  // Declare global and local matrices

  int matrix[M][N];
  int localmat[MP][NP];

  MPI_Comm comm = MPI_COMM_WORLD;

  // Parameters for MPI_Scatterv

  MPI_Datatype vector, resizevector;
  int recvcounts[NPROC], displacements[NPROC];

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  if (size != NPROC)
    {
      if (rank == 0) printf("Must be run on %d processes\n", NPROC);
      MPI_Finalize();

      return;
    }

  if (rank == 0)
    {
      printf("Running on %d processes in %d x %d grid\n", NPROC, ROWP, COLP);

      printf("Global matrix size is %2d x %2d\n", M,  N );
      printf("Local  matrix size is %2d x %2d\n", MP, NP);

      printf("\n");
    }

  ipos = rank / COLP; // row index of this process in the grid
  jpos = rank % COLP; // column index of this process in the grid


  // Initialise local array based on global indices

  for (i=0; i < MP; i++)
    {
      for (j=0; j < NP; j++)
	{
	  iglobal = ipos*MP+i;
	  jglobal = jpos*NP+j;
	    
	  localmat[i][j] = iglobal*N + jglobal;
	}
    }

  // Initialise target global array to -1

  for (i=0; i < M; i++)
    {
      for (j=0; j < N; j++)
	{
	  matrix[i][j] = -1;
	}
    } 

   //  if (rank == 0)
   //    {
   //      printf("matrix on rank 0 before:\n\n");
   //
   //      matprint(&matrix[0][0], M, N);
   //    }
  
   //  printf("Local matrix on rank %d before:\n\n", rank);
   //  matprint(&localmat[0][0], MP, NP);

  MPI_Gather(localmat, MP*NP, MPI_INT, matrix, MP*NP, MPI_INT, 0, comm);

  if (rank == 0)
    {
      printf("matrix on rank 0 after naive gather\n\n");

      matprint(&matrix[0][0], M, N);
      printf("\n");
    }

  // Define datatype for MPxNP subsection of MxN matrix

  MPI_Type_vector(M, NP, N, MPI_INT, &vector);
  MPI_Type_create_resized(vector, 0, sizeof(int), &resizevector);
  MPI_Type_commit(&resizevector);

  // We only receive one datatype from each process

  for (i=0; i < NPROC; i++)
    {
      recvcounts[i] = 1;
    }

  // Compute starting positions of received vectors

  for (i=0; i < ROWP; i++)
    {
      for (j=0; j < COLP; j++)
	{
	  displacements[i*COLP+j] = i*MP*N + j*NP;
	}
    }

  MPI_Gatherv(localmat, MP*NP, MPI_INT,
	      matrix, recvcounts, displacements, resizevector,
	      0, comm);  

  if (rank == 0)
    {
      printf("matrix on rank 0 after vector gather\n\n");

      matprint(&matrix[0][0], M, N);
    }
    
  MPI_Finalize();
}


void matprint(int *mat, int m, int n)
{
  int i, j;

  for (i=0; i < m; i++)
    {
      for (j=0; j < n; j++)
	{
	  printf("%3d ", mat[i*n+j]);
	}
      printf("\n");
    }
}
