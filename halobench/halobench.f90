program benchio

  use benchclock
  use haloswap

  implicit none

  !
  ! Set main paramemers: size of each halo and number of repeitions
  !

  integer, parameter :: nbuf = 10000
  integer, parameter :: nrep = 1000

  integer :: i, ineigh, idim

  double precision, dimension(nbuf, nneigh, ndim) :: sendbuf, recvbuf

  integer :: rank, size, ierr, comm, cartcomm, dblesize
  integer, dimension(ndim) :: dims, coords

  integer, parameter :: mib = 1024*1024

  logical :: reorder = .false.

! Periodic boundaries

  logical, dimension(ndim) :: periods = [.true., .true., .true.]

  double precision :: t0, t1, time, iorate, mibdata

  call MPI_Init(ierr)

  comm = MPI_COMM_WORLD

  call MPI_Comm_size(comm, size, ierr)
  call MPI_Comm_rank(comm, rank, ierr)

  dims(:) = 0

! Set 3D processor grid

  call MPI_Dims_create(size, ndim, dims, ierr)

  call MPI_Type_size(MPI_DOUBLE_PRECISION, dblesize, ierr)

  mibdata = float(dblesize*nbuf*nneigh*ndim)/float(mib)

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Simple haloswap benchmark'
     write(*,*) '-------------------------'
     write(*,*)
     write(*,*) 'Running on ', size, ' process(es)'
     write(*,*) 'Process grid is (', dims(1), ', ', dims(2), ', ', dims(3), ')'
     write(*,*) 'Each halo contains', nbuf, ' doubles'
     write(*,*)
     write(*,*) 'Total amount of halo data = ', mibdata, ' MiB per process'
     write(*,*)
     write(*,*) 'Number of repetitions = ', nrep
     write(*,*) 'Clock resolution is ', benchtick()*1.0e6, ' microseconds'
  end if

  call MPI_Cart_create(comm, ndim, dims, periods, reorder, cartcomm, ierr)

! Set data to illegal values

  recvbuf(:,:,:) = -1
  sendbuf(:,:,:) = -1
  
! Set halo data core to have unique values

  call MPI_Cart_coords(cartcomm, rank, ndim, coords, ierr)
  
  do idim = 1, ndim
     do ineigh = 1, nneigh
        do i = 1, nbuf

           sendbuf(i, ineigh, idim) =   rank*(nbuf*nneigh*ndim) &
                                      + (idim-1)*nbuf*nneigh  &
                                      + (ineigh-1)*nbuf + i

!           write(*,*) 'rank = ', rank, ', x(', i, ', ', ineigh, ', ', idim, &
!                ') = ', sendbuf(i,ineigh,idim)

        end do
     end do
  end do

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Sendrecv'
     write(*,*) '--------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call halosendrecv(nrep, sendbuf, recvbuf, nbuf, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Redblack'
     write(*,*) '--------'
  end if

  if (any(mod(dims(:),2) /= 0)) then

     if (rank == 0) then
        write(*,*) 'Redblack requires all process dimensions to be even'
        write(*,*) 'Skipping this test'
     end if

  else

     call MPI_Barrier(comm, ierr)
     t0 = benchtime()
     call haloredblack(nrep, sendbuf, recvbuf, nbuf, cartcomm)
     call MPI_Barrier(comm, ierr)
     t1 = benchtime()

     time = t1 - t0
     iorate = dble(nrep)*mibdata/time

     if (rank == 0) then
        write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
     end if

  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Isend / Recv / Wait'
     write(*,*) '--------------------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call haloisendrecvwait(nrep, sendbuf, recvbuf, nbuf, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Irecv / Send / Wait'
     write(*,*) '--------------------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call haloirecvsendwait(nrep, sendbuf, recvbuf, nbuf, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Irecv / Isend / Wait (pairwise)'
     write(*,*) '-------------------------------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call haloirecvisendwaitpair(nrep, sendbuf, recvbuf, nbuf, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Irecv / Isend / Waitall'
     write(*,*) '-----------------------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call haloirecvisendwaitall(nrep, sendbuf, recvbuf, nbuf, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Persistent Comms'
     write(*,*) '----------------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call halopersist(nrep, sendbuf, recvbuf, nbuf, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Neighbourhood Collective'
     write(*,*) '------------------------'
  end if

  call MPI_Barrier(comm, ierr)
  t0 = benchtime()
  call haloneighboralltoall(nrep, sendbuf, recvbuf, nbuf, cartcomm)
  call MPI_Barrier(comm, ierr)
  t1 = benchtime()

  time = t1 - t0
  iorate = dble(nrep)*mibdata/time

  if (rank == 0) then
     write(*,*) 'Secs = ', time, ', bwidth = ', iorate, ' MiB/s'
  end if

  if (rank == 0) then
     write(*,*)
     write(*,*) 'Finished'
     write(*,*) '--------'
     write(*,*)
  end if

  call MPI_Finalize(ierr)
  
end program benchio
