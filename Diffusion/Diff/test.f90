   program get_diff

!      
! 
!     This program will calculate diffusivities 
!     from mean square displacement data 
! 
!     author  David Keffer 
!     Department of Chemical Engineering 
!     University of Tennessee, Knoxville 
!     last updated  October 6, 2001 
! 

      implicit double precision (a-h,o-z)   


      integer, parameter :: maxstp = 137300         !number of data production steps 
      integer, parameter :: kmsd = 10               !sampling interval 
      integer, parameter :: N = 64                  !number of molecules 
      double precision, parameter  :: dt = 0.0096756 !timestep (fs) 
      character*12 :: cmsd, cout                    !character variables 
      character*3, dimension(1:4) :: cname 
      double precision, dimension(1:4) :: Dav, Dsd 
      double precision, dimension(1:3) :: slope,slopesd,yinter,yintersd 
      double precision, allocatable :: md_msd(:,:), time_vec(:), & 
        & xmsd(:,:) 
! 
      cout = 'get_diff.out' 
      cmsd = 'h2o-64.mod' 
! number of times represented in data  
      ntime = maxstp/kmsd + 1 
! number of rows of data  
      ndata = N*ntime 
      allocate (md_msd(1:ndata,1:3)) 
      open(unit=1,file=cout,form='formatted',status='unknown') 
      open(unit=2,file=cmsd,form='formatted',status='old') 
      print *, ' ntime = ', ntime, ' ndata = ', ndata 
      do i = 1, ndata, 1 
         read(2,*) md_msd(i,1:3) 
      enddo 
      print *, ' read all the data' 
! number of origins is half number of time steps 
      norigin = (ntime-1)/2 
! minimum number of intervals to contribute to diffusivity 
      nmin = 50 
! maximum number of intervals to contribute to diffusivity 
      nmax = norigin 
! 
      if (nmin .gt. nmax) then 
         print *, ' We have a problem. ' 
         print *, ' nmin = ', nmin, ' nmax = ', nmax 
         stop 
      endif 
! store mean square displacements in xmsd 
      allocate (time_vec(1:norigin), xmsd(1:norigin,1:3)) 
      do i = 1, norigin, 1 
         time_vec(i) = dfloat(i*kmsd)*dt 
      enddo

      xmsd = 0.d0
      do i = 1, N, 1
         do j = 1, norigin, 1
            jstart = (j-1)*N + i
            do k = nmin, nmax, 3
               kend = jstart + k*N
               xmsd(k,1:3) = xmsd(k,1:3) + &
               (md_msd(kend,1:3) - md_msd(jstart,1:3) )**2
            enddo
         enddo
      enddo
      xmsd = xmsd/dfloat(N*norigin)
      do i = 1, 3, 1
         call dllsr(slope(i), slopesd(i), yinter(i), yintersd(i), &
         nmax-nmin+1, time_vec(nmin:nmax), xmsd(nmin:nmax,i))
      enddo
!    report results
      cname(1) = 'x '
      cname(2) = 'y '
      cname(3) = 'z '
      cname(4) = 'avg'
      do i = 1, 3, 1
         write(6,1007) cname(i), slope(i), yinter(i)
         write(1,1007) cname(i), slope(i), yinter(i)
      enddo
1007 format(a3, ' slope = ', e16.8, ' y-intercept = ', e16.8,&
& ' A^2/fs')
      Dav(1:3) = 0.5d0*slope(1:3)*1.0d-5 ! convert to m^2/sec
      Dsd(1:3) = 0.5d0*slopesd(1:3)*1.0d-5 ! convert to m^2/sec
      Dav(4) = sum(Dav(1:3))/3.d0
! standard deviation of average diffusivity
      term1 = 3.d0*(Dav(1)*Dav(1) + Dav(2)*Dav(2) + Dav(3)*Dav(3) )
      Dsd(4) = sqrt( (term1 - Dav(4)*Dav(4)*9.d0) /6.d0 )
      do i = 1, 4, 1
         write(6,1006) cname(i), Dav(i), Dsd(i)
         write(1,1006) cname(i), Dav(i), Dsd(i)
      enddo
1006 format(a3, ' diffusivity avg = ', e16.8, ' stand dev = ', e16.8,&
& ' m^2/sec ')

!
!
!write xmsd vs time data for later plotting
!
!
      do i = 1, norigin, 1
         write(1,1008) time_vec(i), xmsd(i,1:3)
      enddo
1008 format(4(e16.8,1x))
      close (unit=1,status='keep')
      close (unit=2,status='keep')
      stop
      end program

      subroutine dllsr(slope, slopesd, yinter, yintersd, n, x, y)
         implicit double precision (a-h, o-z)
         double precision, intent(out) :: slope, slopesd, yinter, yintersd
         integer, intent(in) :: n
         double precision, intent(in), dimension(1:n) :: x, y
         xn = dfloat(n)
         xavg = sum(x)/xn
         yavg = sum(y)/xn
         sumxy = 0.d0
         sumxx = 0.d0
         sumx2 = 0.d0
         do i = 1, n, 1
            sumxy = sumxy + (x(i) - xavg)*(y(i) - yavg)
            sumxx = sumxx + (x(i) - xavg)*(x(i) - xavg)
            sumx2 = x(i)*x(i)
         enddo
         slope = sumxy/sumxx
         yinter = yavg - slope*xavg
         sse = 0.d0
         do i = 1, n, 1
            sse = sse + (y(i) - slope*x(i) - yinter)**2.d0
         enddo
         sig2 = sse/dfloat(n-2)
         slopesd = dsqrt(sig2/sumxx)
         yintersd = dsqrt(sig2/dfloat(n)*sumx2/sumxx)
         return
      end

