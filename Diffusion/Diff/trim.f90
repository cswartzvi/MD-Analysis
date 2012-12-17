!This Program will take the position output file from QE and  
!trim the unwanted coordinates and the step numbers
!
!
!
!Charles Swartz
!-------------------------------------------------------------

 PROGRAM output_trim

   implicit none

   integer                             :: steps, interval, tot, i, k, &
                                          atoms1, atoms2, step
   real(kind=8), allocatable           :: r(:,:) 
   real                                :: x(3)

   steps    = 137400 
   interval = 10
   !tot is the total number of intervals
   tot      = steps/interval
   atoms1   = 64
   atoms2   = 128

   allocate(r(3,tot))

   open(unit=15, file="h2o-64.pos.mod")

   do i=1,tot,1
      read(*,*) step
      do k=1,atoms1,1
         read(*,*) r(1,i), r(2,i), r(3,i)
         write(15,*) r(1,i), r(2,i), r(3,i)
      enddo
      do k=1,atoms2,1
         read(*,*) x(1), x(2), x(3)
      enddo
   enddo


 END PROGRAM output_trim
