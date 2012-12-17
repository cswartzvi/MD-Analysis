PROGRAM test 
   IMPLICIT NONE

   INTEGER, PARAMETER         :: DP = KIND(0.0d0)  
   INTEGER, PARAMETER         :: maxsp = 20
   REAL, PARAMETER            :: bohr_to_ang = 0.529177249 
   REAL,    DIMENSION(:,:), ALLOCATABLE   :: x(:,:), y(:,:), z(:,:)
   REAL,    DIMENSION(:,:), ALLOCATABLE   :: int_x(:,:), int_y(:,:), int_z(:,:)
   INTEGER, DIMENSION(:), ALLOCATABLE     :: g(:) 
   REAL                    :: xr, yr, zr, r, delg, box, rho, vx, vy, vz, nid,  &
                              dr, rc, bull1, bull2, bullN, vol, xrt,           &
                              yrt, zrt, vb, r2, rc2, box2
   INTEGER                 :: npart, nref1, nref2, npairs, i,j, k, ig, nhis,   &
                              ngr = 0, ios, fcount, fstart, fstop, fstep,      &
                              isave, nsp, nat(maxsp), ref1_start, ref1_end,    &
                              ref2_start, ref2_end, n1, n2, is, read_start,    &
                              read_stop, gplus, ir, ia, label
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
   box = box*bohr_to_ang
   box2 = 2.0*box
   vol = box**3 
   rho = npairs/vol
   nhis = int(rc/dr)
   rc2= rc*rc
   fstart = fstart/isave 
   fstop = fstop/isave
   fstep = fstep/isave
   
   !Alolocate memory
   ALLOCATE(g(nhis))
   ALLOCATE(x(npart*8, nsp), y(npart*8, nsp), z(npart*8, nsp))
   ALLOCATE(int_x(npart, nsp), int_y(npart, nsp), int_z(npart, nsp))


   Do i=1,nhis,1
      g(i)=0
   ENDDO


   !Open and Read in the datafile (starting from sstart->fstop)
   OPEN(15, FILE = TRIM(datafile), STATUS = 'OLD')

   DO fcount=fstart,fstop,fstep

      ngr = ngr + 1
      
      READ(15,*)  bullN, bull1
      read_start = 1
      read_stop  = nat(1)
      DO is=1, nsp
         DO i=read_start, read_stop
            READ(15,*)  int_x(i,is), int_y(i,is), int_z(i,is)
            int_x(i,is) = int_x(i,is)*bohr_to_ang
            int_y(i,is) = int_y(i,is)*bohr_to_ang
            int_z(i,is) = int_z(i,is)*bohr_to_ang
         ENDDO
         read_start  = read_start + nat(is)
         read_stop   = read_stop  + nat(is + 1)
      ENDDO

      !Periodic 
      !DO is=1,nsp,1
      !   DO i=1,nat(is)
      !      int_x(i, is)  = int_x(i, is) - NINT(int_x(i, is)/box)*box
      !      int_y(i, is)  = int_y(i, is) - NINT(int_y(i, is)/box)*box
      !      int_z(i, is)  = int_z(i, is) - NINT(int_z(i, is)/box)*box
      !   END DO
      !END DO

      !The 1st octant: a 2x2x2 unit cell construction for each speices
      !with this construction the algorithm,
      !r = r - INIT(r/(2*box))*2*box 
      !will determine all values in a "box" radius.

      DO is =1,nsp,1
         !This with be the "true" atomic counter, for the first octant
         ir = 0
         DO i=0,1,1
            DO j=0,1,1
               DO k=0,1,1
                  DO ia=1,nat(is)
                     ir = ir + 1
                     x(ir, is) = int_x(ia, is) + i*box
                     y(ir, is) = int_y(ia, is) + j*box
                     z(ir, is) = int_z(ia, is) + k*box
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      !Determine if nref1 == nref2
      !Which "type" of calculation is to be preformed
      if(nref1 == nref2) then
         ref1_end    =  nat(nref1)*8 - 1     
         gplus       =  2
      else
         ref1_end    =  nat(nref1)*8
         gplus       =  1
      endif

      DO n1=1,ref1_end,1

         if(nref1 == nref2) then
            ref2_start  =  n1 + 1
         else
            ref2_start  =  1
         endif

         DO n2=ref2_start, (nat(nref2)*8), 1


                  xr = abs(x(n1, nref1) - x(n2, nref2))
                  yr = abs(y(n1, nref1) - y(n2, nref2))
                  zr = abs(z(n1, nref1) - z(n2, nref2))

                  xr = xr - NINT(xr/(box2))*box2
                  yr = yr - NINT(yr/(box2))*box2
                  zr = zr - NINT(zr/(box2))*box2

                  r2 = xr**2 + yr**2 + zr**2

                  if(r2 < (rc2)) then
                  ig = nint(sqrt(r2)/dr)
                     g(ig) = g(ig) + gplus
                  endif


         ENDDO
      ENDDO

   ENDDO

   CLOSE(15)

   do i=1,nhis,1
      r=dr*(i-1)
      vb=(i**3-(i-1)**3)*dr**3
      nid= (4./3.)*pi*vb
         
      print *, r, real(g(i)/(8*ngr*rho*nid)) 
   enddo

   DEALLOCATE(x, y, z, g)

END PROGRAM test
