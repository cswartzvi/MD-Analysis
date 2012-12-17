 Program inverse

   implicit none

   real        :: a(3,3), b(3,3), c(3,3), den
   integer     :: i, j

   a(1,1:3) = (/1, 2, 1/)
   a(2,1:3) = (/3, 5, 6/)
   a(3,1:3) = (/7, 8, 9/)

   call invert(a, b, den)
   
   do i=1,3,1
         write(*,1000) (a(i,j),j=1,3)
   enddo
   
   write(*,*) ''
   write(*,*) den
   write(*,*) ''

   do i=1,3,1
         write(*,1000) (b(i,j),j=1,3)
   enddo

   c = MATMUL(a, b)

   print *, ''

   do i=1,3,1
         write(*,1000) (c(i,j),j=1,3)
   enddo

 1000 Format ( 3(F10.5) )
   

   contains

      subroutine invert(Mi, Mo, det)

         implicit none

         real, intent(in)     :: Mi(3,3)
         real, intent(out)    :: Mo(3,3), det
         real                 :: tmp(3,3), s

         integer              :: i,j,k,l  !indexes
         integer              :: n,ir     !int counters

         det = 0.0
         s   = 1.0
         i   = 1
         j   = 2
         k   = 3


         do
            do n=1,3,1
               det= det + s*Mi(1,i)*Mi(2,j)*Mi(3,k)
               l = i
               i = j
               j = k
               k = l
            end do

            i = 2
            j = 1
            k = 3
            s = -s
            if (s .GT. 0.0) then
               exit     
            endif
         enddo

         IF(ABS(det) .LT. 1.0e-20) then
            write(*,*) 'Error: Singular Matrix'
            Stop   
         endif

         i = 1
         j = 2
         k = 3

         do ir=1,3
            tmp(ir,1) = (Mi(2,j)*Mi(3,k) - Mi(2,k)*Mi(3,j)) / det
            tmp(ir,2) = (Mi(3,j)*Mi(1,k) - Mi(3,k)*Mi(1,j)) / det
            tmp(ir,3) = (Mi(1,j)*Mi(2,k) - Mi(1,k)*Mi(2,j)) / det
         
            l = i
            i = j
            j = k
            k = l
         enddo
        
         do l=1,3
            do k=1,3
               Mo(k,l) = tmp(k,l)
            enddo
         enddo 

      end subroutine invert

 End Program inverse
