module trafficlib

  use unirand

  implicit none

  integer, parameter :: seed = 5743

contains

integer function initroad(road, n, density, seedval)

  integer :: n, seedval
  integer :: road(n)
  real    :: density
  
  integer :: i, ncar
  real    :: rng

  call rinit(seedval)

  ncar = 0

  do i = 1, n

     rng = uni()

     if (rng .lt. density) then

        road(i) = 1

     else

        road(i) = 0

     end if

     ncar = ncar + road(i)
    
  end do

  initroad = ncar

end function initroad


integer function updateroad(newroad, oldroad, n)

  integer :: n
  integer :: newroad(0:n+1), oldroad(0:n+1)

  integer :: i, nmove

  nmove = 0

  do i = 1, n

     if (oldroad(i) .eq. 1) then

        if (oldroad(i+1) .eq. 1) then

           newroad(i) = 1
           
        else

           newroad(i) = 0
           nmove = nmove + 1

        end if

     else

        if (oldroad(i-1) .eq. 1) then

           newroad(i) = 1

        else

           newroad(i) = 0
           
        end if

     end if

  end do

  updateroad = nmove

end function updateroad

end module trafficlib

