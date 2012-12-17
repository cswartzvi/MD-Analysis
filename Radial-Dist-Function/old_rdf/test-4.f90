PROGRAM test 
   IMPLICIT NONE

   INTEGER, PARAMETER         :: DP = KIND(0.0d0)  
   INTEGER, PARAMETER         :: maxsp = 20
   REAL, PARAMETER            :: bohr_to_ang = 0.529177249 
   REAL,    DIMENSION(:), ALLOCATABLE     :: x(:), y(:), z(:)
   INTEGER, DIMENSION(:), ALLOCATABLE     :: g(:) 
   REAL                    :: xr, yr, zr, r, delg, box, rho, vx, vy, vz, nid,  &
                              dr, rc, bull1, bull2, bullN, vol, xrt,      &
                              yrt, zrt, vb, r2
   INTEGER                 :: npart, nref1, nref2, npairs, i,j, k, ig, nhis,   &
                              ngr = 0, ios, fcount, fstart, fstop, fstep,      &
                              isave, nsp, nat(maxsp), ref1_start, ref1_end,    &
                              ref2_start, ref2_end, n1, n2
   REAL, PARAMETER         :: pi = 3.14159
   CHARACTER(len=200)   :: datafile
   CHARACTER(len=4)     :: format_string = ".xyz"
   CHARACTER(len=500)   :: x1 
   LOGICAL              :: first = .TRUE., eqpairs = .FALSE.

   !Initial Values
   NAMELIST /input_rdf/ datafile, fstart, fstop, fstep, isave, npart, box,     &
      rc, dr, nsp, nat, nref1, nref2
   READ(*, input_rdf, iostat=ios)

   !Initial Calculations
   if (nref1 == nref2) THEN
      npairs = nat(nref1)*(nat(nref2)-1)
      eqpairs=.TRUE.
   else
      npairs = nat(nref1)*nat(nref2)
   endif
   box = box*(0.529177249)
   vol = box**3 
   rho = npairs/vol
   nhis = int(rc/dr)   
   fstart = fstart/isave 
   fstop = fstop/isave
   fstep = fstep/isave
   
   !Alolocate memory
   ALLOCATE(x(npart), y(npart), z(npart), g(nhis))

   DO i=1,npart,1
      x(i)=0.0
      y(i)=0.0
      z(i)=0.0
   ENDDO

   Do i=1,nhis,1
      g(i)=0
   ENDDO


   !Determine particle reference labels
   ref1_start = 1
   ref2_start = 1
   ref1_end   = nat(1)
   ref2_end   = nat(1)

   IF(nref1 > 1) THEN
      DO i=2,nref1,1
         ref1_start = ref1_start + nat(i-1)
      ENDDO
   ENDIF
   IF(nref2 > 1) THEN
      DO i=2,nref2,1
         ref2_start = ref2_start + nat(i-1)
      ENDDO
   ENDIF

   !Open and Read in the datafile (starting from sstart->fstop)
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
  
 

      !Detetmine if the piars of ref1 and ref2 are the same 
!      if (testing) then
        do n1=ref1_start,(ref1_end-1),1
            do n2=(n1+1),ref2_end,1
!         do n1=ref1_start, ref1_end, 1
!            do n2=ref2_start, ref2_end, 1

               if(n1 /= n2) then
               xr = abs(x(n1) - x(n2))
               !xrt = xr - box
               xr = xr - box*nint(xr/box)

               yr = abs(y(n1) - y(n2))
               !yrt = yr - box 
               yr = yr - box*nint(yr/box)

               zr = abs(z(n1) - z(n2))
               !zrt = zr - box 
               zr = zr - box*nint(zr/box)
            

               r2 = xr**2 + yr**2 + zr**2

               if(r2 < (rc*rc)) then
               ig = nint(sqrt(r2)/dr)
                  g(ig) = g(ig) +2 
               endif
   
               !r = sqrt(xrt**2 + yrt**2 + zrt**2)

               !ig = nint(r/dr)
               !if(r < rc) then
               !   g(ig) = g(ig) + 2
               !endif
               !2 Due to ith AND jth particle!!
               endif

               


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
