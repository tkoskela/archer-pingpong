program mpigather2d

  !
  ! This extends the 1D "cyclicgather" example to a 2D distributed array.
  ! For example, consider a 4x4 matrix distribute across a 2x2 process grid:
  !
  !     Matrix             Process grid
  !
  !  0  4 |  8 12      
  !  1  5 |  0 13        rank 0 | rank 1
  !  -----+------        -------+-------
  !  2  6 | 10 14        rank 2 | rank 3
  !  3  7 | 11 15 
  !
  ! Use a combination of MPI_Gatherv and a resized vector datatype to
  ! reassemble this in the correct order on rank 0.
  !
  ! Note we use C-based ordering for process coordinates in the grid
  ! as this is the convention used by MPI in Cartesian topologies.
  !
  
  use mpi
  
  implicit none
  
  integer, parameter :: M=12    ! Global number of rows
  integer, parameter :: N=10    ! Global number of columns

  integer, parameter :: ROWP=3  ! Number of processes for row distribution
  integer, parameter :: COLP=2  ! Number of processes for column distribution

  integer, parameter :: NPROC=ROWP*COLP  ! Total number of processes

  integer, parameter :: MP=M/ROWP  ! Local number of rows
  integer, parameter :: NP=N/COLP  ! Local number of columns

  integer :: i, j, rank, size, ierr
  integer :: iglobal, jglobal, ipos, jpos

  ! Declare global and local matrices

  integer, dimension(M, N)  :: matrix
  integer, dimension(MP,NP) :: localmat

  integer :: comm = MPI_COMM_WORLD

  call MPI_Init(ierr)
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, size, ierr)

  if (size /= NPROC) then

     if (rank == 0) write(*,*) 'Must be run on ', NPROC, ' processes'
     call MPI_Finalize(ierr)

     stop

  end if

  if (rank == 0) then

     write(*,*) 'Running on ', NPROC, ' processes in ', &
                 ROWP, ' x ', COLP, ' grid'
     write(*,*) 'Global matrix size is ', M,  ' x ', N
     write(*,*) 'Local  matrix size is ', MP, ' x ', NP
     write(*,*)

  end if

  ipos = rank / COLP      ! row index of this process in the grid
  jpos = mod(rank, COLP)  ! column index of this process in the grid


  ! Initialise local array based on global indices

  do i = 1, MP
     do j = 1, NP

        iglobal = ipos*MP+i;
        jglobal = jpos*NP+j;
	    
        localmat(i,j) = (jglobal-1)*M + iglobal;

     end do
  end do

  ! Initialise target global array to -1

  do i = 1, M
     do j = 1, N

        matrix(i,j) = -1

     end do
  end do

  call MPI_Gather(localmat, MP*NP, MPI_INTEGER, &
                  matrix, MP*NP, MPI_INTEGER, 0, comm, ierr)

  if (rank == 0) then

     write(*,*) 'matrix on rank 0 after naive gather'
     write(*,*)

     call matprint(matrix, M, N);
     write(*,*)
     
  end if

  call MPI_Finalize(ierr)

end program mpigather2d


subroutine matprint(mat, m, n)

  implicit none

  integer :: m, n
  integer, dimension(m,n) :: mat

  integer :: i, j
  
  do i = 1, m
     write(*,fmt='(99(i3,1x))') (mat(i,j), j = 1, N)
  end do

end subroutine matprint
