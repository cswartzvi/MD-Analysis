 Program Diff

   implicit none


   integer,parameter    :: spmax=5

   integer              :: maxstep,    &  !Max Number of production steps (actual)
                           nsp(spmax), &  !Total number of species, for each type nspt
                           mint,       &  !Sampling interval
                           nspt,       &  !Total number of species
                           aindex,     &  !Atom index to be considered displacement, tagged
                           Nat,        &  !Total number of atoms for aindex
                           nsample,    &  !Total Number of samples
                           ndata,      &  !Total Number of datapoints
                           norigin,    &  !Total number of origins to consider in the average
                           nmin,nmax,  &  !Min/Max intervals for each origin's average
                           i,j,k,n,    &  !simple indexes   
                           jstart,     &
                           kend
                                

   real, allocatable    :: tau(:,:),   &  !Atom positions tau(atom-index,dim), for aindex
                                          ! Irrelevant species will be trimed as
                                          ! *.pos file is read 
                           msd(:),     &  !Total mean squared displacement
                           time(:)        !time vector

   real                 :: dx, dy, dz, &  !Displacements 
                           dt,         &  !time step for ONE interval (in fs)
                           dummy,      &  !A dummy variable
                           dum(3)         !A dummy array


   character(len=15)     :: posfile, outfile

   

   Namelist /input/ maxstep, mint, dt, nspt, nsp, aindex, posfile, outfile

   read(*,input)

   !Initial Values 
   maxstep = maxstep - mint 
   nsample = maxstep/mint + 1
   ndata =  nsp(aindex)*nsample
   
   call print_out(0)

   !Number of Origins is half of the total nummber of samples (configurations)
   !norigin = (nsample-1)/2
   norigin = (nsample-1)/4
   !Number of min max samples
   nmin  = 50 
   !nmin = 1
   !nmax  =  norigin
   nmax = 3*(nsample-1)/4
   !Number of atoms in consideration (tagged)
   Nat = nsp(aindex)

   if(nmin > nmax) then
      write(*,*) ' ERROR: Not enough origins! '
      write(*,*)  ' nmin = ', nmin,  'nmax = ', nmax
      stop
   endif

   !allocate( tau(1:ndata,1:3), time(1:norigin), msd(1:norigin) )
   allocate( tau(1:ndata,1:3), time(1:nmax), msd(1:nmax) )

   !Open Files
   open(unit=1,file=outfile,form='formatted',status='unknown')
   open(unit=2,file=posfile,form='formatted',status='old')

   !Read in the correct the atoms aindex
   n = 0
   do i=1,nsample,1
      read(2,*) dummy, dummy
      do j=1,nspt,1
         do k=1,nsp(j),1
            if(j==aindex) then
               n = n + 1
               read(2,*) tau(n,1:3)
               write(99,*) tau(n,1:3)
            else
               read(2,*) dum(1:3)
            endif 
         enddo
      enddo
   enddo
   call print_out(1)

   !calculate the time vector in picosecond
   do i=1,nmax,1
      time(i) = real(i*mint)*dt*0.001
   enddo
   
   msd(1:norigin) = 0.0
   
   call print_out(2)
   do i =1,Nat,1
      !Loop over the Origins (Atoms in different configurations)
      do j=1,norigin,1
         !Set the Origin to the "Next" atom
         jstart = (j-1)*Nat + i
         !Loop over the intervals
         do k=nmin,nmax,1
            kend = jstart + k*Nat
            dx = tau(kend,1) - tau(jstart,1)
            dy = tau(kend,2) - tau(jstart,2)
            dz = tau(kend,3) - tau(jstart,3)
            msd(k) = msd(k) + dx**2 + dy**2 + dz**2
         enddo   
      enddo
   enddo
   call print_out(3)

   msd(:) = msd/(real(Nat*norigin))

   call print_out(11)   
 
   call print_out(4)

   deallocate(time, msd, tau)

   contains
      
      !---------------------------------------------------------------------------------------
      !Print-out Subroutine
      !io selects from a variety of output functions
      !---------------------------------------------------------------------------------------
      Subroutine print_out(io)
   
         implicit none

         integer, intent(in)  :: io

         select case (io)
            case(0)
                  write(*,*)
                  write(*,*)  '---------------------------------------------------------'
                  write(*,*)  '          Mean Sqaured Displacment Program               '
                  write(*,*)  '               Charles W. Swartz VI                      '
                  write(*,*)  '---------------------------------------------------------'
                  write(*,*)  ' nsample = ', nsample, ' ndata = ', ndata
            case(1)
                  write(*,*)  ' Data from ', trim(posfile), ' read in correctly.'
            case(2)
                  write(*,*)
                  write(*,*) ' Begining main loop (this may take a few moments) ...'
            case(3)
                  write(*,*) ' ... Main loop complete!'
            case(4)
                  write(*,*)
                  write(*,*) ' MSD data output-file : ', outfile
                  write(*,*) ' Program Complete!'
                  write(*,*)
            case(11)
                  write(1,*) '#This file contains time and MSD values for ',posfile
                  do i=1,nmax,1
                     write(1,100) time(i), msd(i)*(0.52917)**2
                  enddo

            case default 
         end select

100 FORMAT (2(e14.8,2x))

      end subroutine print_out

 End Program Diff
