!
!
!-----------------------------------------------------------------------------------------
PROGRAM cp_dos
   !--------------------------------------------------------------------------------------
   !Notes:
   !  This program will perform the calculations with enegery in RY, and
   !  AFTERWARDS switch back to eV
   !
   !
   !
   !--------------------------------------------------------------------------------------
   implicit none
   !
   integer, parameter   :: DP  = SELECTED_REAL_KIND(15,99)
   integer, parameter   :: ncol              = 10
   REAL(DP), PARAMETER  :: pi                = 3.14159265358979323846_DP
   REAL(DP), PARAMETER  :: sqrtpi            = 1.77245385090551602729_DP
   REAL(DP), PARAMETER  :: sqrtpm1           = 1.0_DP / sqrtpi
   REAL(DP), PARAMETER  :: HARTREE_SI        = 4.35974394E-18_DP   ! J
   REAL(DP), PARAMETER  :: ELECTRONVOLT_SI   = 1.602176487E-19_DP  ! J
   REAL(DP), PARAMETER  :: AUTOEV            = HARTREE_SI / ELECTRONVOLT_SI
   REAL(DP), PARAMETER  :: rytoev            = AUTOEV / 2.0_DP
   !------------------------------------------------------
   !
   integer                 :: step, stepNum, nbnd, ndos, ios, nrow, remain, i, j, n
   
   real(DP)                :: E, Emin, Emax, dos, DeltaE, degauss, time, wk

   real(DP), allocatable   :: eig(:)

   character(len=256)      :: dosfile, eigfile
   character(len=256)      :: dummy 

   NAMELIST /indos/ dosfile, eigfile, degauss, DeltaE, step, nbnd, wk, Emin, Emax

   READ (*, nml=indos, iostat=ios )
   if (ios /=  0) then
   write(*,*) ""
      write(*,*) " Required values:"
      write(*,*) "    eigfile:   Previous CP calculation eigenvalues (*.eig file)"
      write(*,*) "    filedos:   Desired output file"
      write(*,*) "    step:      Configuration snapshot number"
      write(*,*) "    DeltaE:    Energy grid step (eV)"
      write(*,*) "    Emin:      Min. Energy Scale (eV)"
      write(*,*) "    Emax:      Max Energy Scale (eV)"
      write(*,*) "    degauss:   Gaussian broadening, Ry (not eV!)"
      write(*,*) "    nbnd:      Total number of bands" 
      write(*,*) "    wt:        Weight for the Gamma kpoint"
      write(*,*) " -------------------------"
      write(*,*) ""
      stop
   endif
   
   nrow = CEILING(nbnd/real(ncol))
   remain = ncol - (nrow*ncol - nbnd)
   !write(*,*) ' nrow: ', nrow
   !write(*,*) ' remain: ', remain
   allocate ( eig(nbnd) )

   !write(*,*) nrow, ncol, remain, nbnd
   
   !open eigfile
   open(unit=1, file=TRIM(eigfile), iostat=ios)
   if (ios /= 0) then
      write(*,*) " ERROR reading eigfile: ", TRIM(eigfile), ". Exiting program..."
      stop
   endif
   !open dosfile
   open(unit=2, file=TRIM(dosfile), iostat=ios)
   if (ios /= 0) then
      write(*,*) " ERROR reading dosfile: ", TRIM(dosfile), ". Exiting program..."
      stop
   endif

   !Read In the *.eig
   readloop: do
   read(1,*) dummy, stepNum, time
   !write(*,*) TRIM(dummy), stepNum, time
      if( stepNum == step) then
         write(*,*) '    Configuration Found, step : ', step
         read(1,'(A)') dummy
         do i=1,(nrow-1)
            read(1,*) (eig((i-1)*ncol + j), j=1,ncol)
         enddo
         read(1,*) (eig((nrow-1)*ncol + j), j=1,remain)
         exit
      else
         do i=1,(nrow+1)
            read(1,'(A)') dummy
         enddo
      endif
   enddo readloop
 
   !Record the read in eignevalues in fort.4
   do i=1,nbnd
      write(4,'(f8.3)') eig(i)
   enddo

   !Emax = eig(nbnd)/rytoev + 3.0*degauss
   Emax = Emax/rytoev 
   Emin = Emin/rytoev
   DeltaE = DeltaE/rytoev

   ndos = NINT ( (Emax - Emin) / DeltaE+0.500001d0)

  
   !WRITE(2,'("#  E (eV)      dos(E) ")')
   do n=1, ndos
      dos = 0.0d0
      E = Emin + (n - 1) * DeltaE

      do i= 1, nbnd
         dos  = dos  + wk  * w0gauss ( (E- eig(i)/rytoev)/ degauss)
         !write(*,*) dos
      enddo

      dos = dos / degauss

      WRITE (2, '(2x, f7.3,3e12.4)') E * rytoev, dos/rytoev
   enddo

   close(1)
   close(2)

 CONTAINS

   function w0gauss (x)

         !delta function

         implicit none
         real(DP) :: w0gauss, x
         ! output: the value of the function
         ! input: the point where to compute the function

         real(DP) :: arg

         arg = min (200.d0, x**2)
         w0gauss = exp ( - arg) * sqrtpm1
      
         return
   end function w0gauss

END PROGRAM cp_dos
