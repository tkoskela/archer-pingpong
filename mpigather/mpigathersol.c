#include <stdio.h>
#include <mpi.h>

/*
 * This exercise starts of with data distributed in a cyclic manner
 * across processes.
 *
 * For example, for NBUFF = 4 running on 3 processes, the sendbuff
 * array is distributed as:
 *
 *     rank 0         rank 1       rank 2
 *
 *   0  3  6  9     1  4  7 10    2  5  8 11
 *
 * We want to use MPI_Gather to bring this back together on to recvbuff
 * on rank 0. We *want* the data to be in the natural order:
 *
 *   recvbuff:  0  1  2  3  4  5  6  7  8  9 10 11
 *
 * Without using MPI_Datatypes, however, the straightforward MPI_Gather
 * call gives data in this order:
 *
 *   recvbuff:  0  3  6  9  1  4  7 10  2  5  8 11
 *
 * The exercise is to use a vector datatype at the receive side to place
 * the data directly into the correct locations in recvbuff. This requires
 * the vector type to be both defined *and resized* appropriately.
 */

#define NBUFF 4

void main(void)
{
  int i;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Datatype vector, resizevector;
  int rank, size;
  
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);

  int sendbuff[NBUFF];
  int recvbuff[size*NBUFF];

  for (i=0; i < NBUFF; i++)
    {
      sendbuff[i] = rank+i*size;
    }

  for (i=0; i < NBUFF*size; i++)
    {
      recvbuff[i] = -1;
    }

  MPI_Gather(sendbuff, NBUFF, MPI_INT, recvbuff, NBUFF, MPI_INT, 0, comm);

  if (rank == 0)
    {
      printf("Result with standard gather\n\n");

      for(i=0; i < NBUFF*size; i++)
	{
	  printf("%d ", recvbuff[i]);
	}
      printf("\n");
    }

  MPI_Type_vector(NBUFF, 1, size, MPI_INT, &vector);
  MPI_Type_create_resized(vector, 0, sizeof(int), &resizevector);
  MPI_Type_commit(&resizevector);


  MPI_Gather(sendbuff, NBUFF, MPI_INT, recvbuff, 1, resizevector, 0, comm);

  if (rank == 0)
    {
      printf("\nResult with vector gather\n\n");

      for(i=0; i < NBUFF*size; i++)
	{
	  printf("%d ", recvbuff[i]);
	}
      printf("\n");
    }

  MPI_Finalize();
}
