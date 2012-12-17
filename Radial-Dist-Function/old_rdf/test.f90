PROGRAM test 
 
   IMPLICIT NONE

   REAL,    DIMENSION(:), ALLOCATABLE     :: x(:), y(:), z(:)
   INTEGER, DIMENSION(:), ALLOCATABLE     :: g(:)  
   REAL                    :: xr, yr, zr, r, delg, box, rho, vx, vy, vz, nid, vb
   REAL                    :: dr, rc, rc2
   INTEGER                 :: npart, bull1, bull2, bullN, i, ig, j, nhis, ngr = 0
   INTEGER                 :: fcount, fstart, fstop, fstep
   REAL, PARAMETER         :: pi = 3.14159
   CHARACTER(len=200)   :: filename
   CHARACTER(len=4)     :: format_string = ".xyz"
   CHARACTER(len=500)   :: x1 
   LOGICAL              :: first=.TRUE.

   fstart = 10000
   fstop  = 590000
   fstep = 1000
 
   npart = 108
   rho = 0.85
   rc = 2.501
   dr = 0.02 

   
   rc2 = rc*rc
   box = (npart/rho)**(0.333333333)
   nhis = int(rc/dr)   

   ALLOCATE(x(npart))
   ALLOCATE(y(npart))
   ALLOCATE(z(npart))

   DO i=1,npart,1
      x(i)=0.0
      y(i)=0.0
      z(i)=0.0
   ENDDO

   ALLOCATE(g(nhis))

   Do i=1,nhis,1
      g(i)=0
   ENDDO

   DO fcount=fstart,fstop,fstep

      ngr = ngr + 1
      
      WRITE(x1,"(I0)") fcount  !Nice example of how to turn a integer into a
                               !string with 'internal files!!!! 

      filename = TRIM(x1)//format_string
      WRITE(13,*) "../../Classical_MD/"//TRIM(filename)

      OPEN(15, FILE = "../../Classical_MD/"//TRIM(filename), STATUS = 'OLD')
      

      READ(15,*)  bullN, bull1
      DO i=1, npart
          READ(15,*) bull2, x(i), y(i), z(i), vx, vy, vz
      ENDDO
   
      CLOSE(15)

      do i=1,(npart-1),1
         do j=(i+1),npart,1

            xr = abs(x(i) - x(j))
            xr = box*nint(xr/box) - xr
            yr = abs(y(i) - y(j))
            yr = box*nint(yr/box) - yr
            zr = abs(z(i) - z(j))
            zr = box*nint(zr/box) - zr

            r = sqrt(xr**2 + yr**2 + zr**2)

            ig = nint(r/dr)
            if(r < rc) then
               g(ig) = g(ig) + 2
            endif
            !2 Due to ith AND jth particle!!


         enddo
      enddo
   
   ENDDO
   

   do i=1,nhis,1
      r=dr*(i-1)
      vb=(i**3-(i-1)**3)*dr**3
      nid= (4./3.)*pi*vb*rho
         
      print *, r, real(g(i)/(ngr*(npart-1)*nid)) 
   enddo

   DEALLOCATE(x, y, z, g)

END PROGRAM test
