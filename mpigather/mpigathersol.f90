program mpigather

  !
  ! This exercise starts of with data distributed in a cyclic manner
  ! across processes.
  !
  ! For example, for NBUFF = 4 running on 3 processes, the sendbuff
  ! array is distributed as:
  !
  !     rank 0         rank 1       rank 2
  !
  !   0  3  6  9     1  4  7 10    2  5  8 11
  !
  ! We want to use MPI_Gather to bring this back together on to recvbuff
  ! on rank 0. We *want* the data to be in the natural order:
  !
  !   recvbuff:  0  1  2  3  4  5  6  7  8  9 10 11
  !
  ! Without using MPI_Datatypes, however, the straightforward MPI_Gather
  ! call gives data in this order:
  !
  !   recvbuff:  0  3  6  9  1  4  7 10  2  5  8 11
  !
  ! The exercise is to use a vector datatype at the receive side to place
  ! the data directly into the correct locations in recvbuff. This requires
  ! the vector type to be both defined *and resized* appropriately.
  !
  
  use mpi
  
  implicit none
  
  integer, parameter :: NBUFF = 4

  integer :: i, rank, size, comm, ierr

  integer :: vectortype, resizevectortype
  integer :: intsize
  integer(kind=MPI_ADDRESS_KIND) :: newsize, offset
  
  integer :: sendbuff(NBUFF)
  integer, allocatable :: recvbuff(:)

  comm = MPI_COMM_WORLD

  call MPI_Init(ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, size, ierr)

  allocate(recvbuff(size*NBUFF))

  do i = 1, NBUFF
     sendbuff(i) = rank + (i-1)*size
  end do

  recvbuff(:) = -1

  call MPI_Gather(sendbuff, NBUFF, MPI_INT, &
                  recvbuff, NBUFF, MPI_INT, &
                  0, comm, ierr)

  if (rank == 0) then

     write(*,*) 'Result with standard gather'
     write(*,*)
     write(*,fmt='(99(i3,1x))') (recvbuff(i), i = 1, NBUFF*size)

  end if

  call MPI_Type_size(MPI_INTEGER, intsize, ierr)

  offset  = 0
  newsize = intsize

  call MPI_Type_vector(NBUFF, 1, size, MPI_INTEGER, vectortype, ierr)

  call MPI_Type_create_resized(vectortype, offset, newsize, &
                               resizevectortype, ierr)

  call MPI_Type_commit(resizevectortype, ierr)

  call MPI_Gather(sendbuff, NBUFF, MPI_INTEGER,  &
                  recvbuff, 1, resizevectortype, &
                  0, comm, ierr);

  if (rank == 0) then

     write(*,*)
     write(*,*) 'Result with vector gather'
     write(*,*)
     write(*,fmt='(99(i3,1x))') (recvbuff(i), i = 1, NBUFF*size)

  end if

  call MPI_Finalize(ierr)

end program mpigather
