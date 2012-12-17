 PROGRAM evp_reader

   implicit none

   integer     :: i, j, tot = 493
   real        :: value(10), tempTot=0.0, avg
   
   do j=1,10,1
      value(j) = 0.0
   enddo

   do i=1,tot,1
      read(*,*) value(:)
      tempTot = tempTot + value(4) 
   enddo

   print *, 'tempTot = ', tempTot
   
   avg = tempTot/tot

   write(*,*) 'Average temperature =', avg



 END PROGRAM evp_reader
