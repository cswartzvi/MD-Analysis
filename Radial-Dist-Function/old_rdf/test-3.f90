PROGRAM test 
   IMPLICIT NONE

   INTEGER, PARAMETER         :: DP = KIND(0.0d0)   
   REAL,    DIMENSION(:), ALLOCATABLE     :: x(:), y(:), z(:)
   INTEGER, DIMENSION(:), ALLOCATABLE     :: g(:)
   REAL  
   REAL                    :: xr, yr, zr, r, delg, box, rho, vx, vy, vz, nid, vb
   REAL                    :: dr, rc, rc2, bull1, bull2, bullN, vol, xrt, yrt
   REAL                    :: zrt
   INTEGER                 :: npart, nref1, nref2, npairs, i, ig, j, nhis, ngr = 0
   INTEGER                 :: ios, fcount, fstart, fstop, fstep, isave, nsp
   REAL, PARAMETER         :: pi = 3.14159
   CHARACTER(len=200)   :: datafile
   CHARACTER(len=4)     :: format_string = ".xyz"
   CHARACTER(len=500)   :: x1 
   LOGICAL              :: first = .TRUE., testing = .TRUE.

   !Initial Values
   NAMELIST /input_rdf/ datafile, fstart, fstop, fstep, isave, npart,  box, &
      rc, dr, nsp
   READ(*, input_rdf, iostat=ios)
   !fstart = 200 
   !fstop  = 65600
   !fstep = 50
   !isave = 50
   !npart = 96
   !nref1 = 32
   !nref2 = 32
   !box = 18.66558 
   !rc = 8.0
   !dr = 0.02 
   !print *, TRIM(datafile)

   !Initial Calculation
   box = box*(0.529177249)
   vol = box**3 
   npairs = nref1*(nref1-1)
   rho = npairs/vol
   nhis = int(rc/dr)   
   fstart = fstart/isave 
   fstop = fstop/isave
   fstep = fstep/isave

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

   OPEN(15, FILE = TRIM(datafile), STATUS = 'OLD')

   DO fcount=fstart,fstop,fstep

      ngr = ngr + 1
      
      READ(15,*)  bullN, bull1
      DO i=1, npart
         READ(15,*)  x(i), y(i), z(i)
         x(i) = x(i)*(0.529177249)
         y(i) = y(i)*(0.529177249)
         z(i) = z(i)*(0.529177249)
      ENDDO
  
!CREATE; ref start and ref end based on the number of nsp in ref1 and ref1
 

      !Detetmine if the piars of ref1 and ref2 are the same 
!      if (testing) then
         do i=1,(nref1-1),1
            do j=(i+1),nref1,1

               xr = abs(x(i) - x(j))
               xrt = xr - box
               yr = abs(y(i) - y(j))
               yrt = yr - box 
               zr = abs(z(i) - z(j))
               zrt = zr - box 
            

               r = sqrt(xr**2 + yr**2 + zr**2)

               ig = nint(r/dr)
               if(r < rc) then
                  g(ig) = g(ig) + 2 
               endif
   
               r = sqrt(xrt**2 + yrt**2 + zrt**2)

               ig = nint(r/dr)
               if(r < rc) then
                  g(ig) = g(ig) + 2
               endif
               !2 Due to ith AND jth particle!!


               


            enddo
         enddo

   ENDDO

   CLOSE(15)

   do i=1,nhis,1
      r=dr*(i-1)
      vb=(i**3-(i-1)**3)*dr**3
      nid= (4./3.)*pi*vb
         
      print *, r, real(g(i)/(ngr*rho*nid)) 
   enddo

   DEALLOCATE(x, y, z, g)

END PROGRAM test
