 program diff
!=======================================================================
!===    Calcola il MSD e il displacement per particella.
!===    Ha bisogno in input del file di traiettorie continue degli atomi  
!===    in coordinate reali.
!===    Output files: fort.2   disp. max per particella specie a 
!===                  fort.12  disp. max per particella specie b 
!===                  fort.3   disp. per particella specie a
!===                  fort.13  disp. per particella specie b
!===                  fort.4   MSD   specie a , b
!===                  
!===    Ultima modifica Carlo Cavazzoni. 
!=======================================================================

      IMPLICIT NONE

      real*8, parameter :: picosecond = 0.2418901D-4


      real*8, allocatable :: taui(:,:,:)
      real*8, allocatable :: tau(:,:,:)
      real*8, allocatable :: sqd(:,:)
      real*8, allocatable :: msqd(:)
      real*8              :: cdm(3),cdmi(3), r(3)
      real*8, allocatable :: mass(:)
      real*8              :: dt, masst, toffset, time
      integer, allocatable :: na(:)
      integer              :: nsp, nstep, nskip, nat, nax
      integer              :: is,ia,k,t,ierr
      character*80         :: atofile, dummy

!=======================================================================

      read(5,*,IOSTAT=ierr) atofile,nstep,nskip,dt,nsp
      if(ierr.ne.0) then
        write(6,*) 'usage: msd.x < input'
        write(6,*) '  input layout'
        write(6,*) '  file.pos nstep nskip timestep nsp'
        write(6,*) '  na(1) mass(1)'
        write(6,*) '  .'
        write(6,*) '  .'
        write(6,*) '  na(nsp) mass(nsp)'
        stop
      endif
 
      toffset = 0.0d0
      allocate(na(nsp))
      allocate(mass(nsp))

      do is = 1, nsp
        read(5,*) na(is),mass(is)
      end do

!=======================================================================

      nat   = sum(na)
      masst = sum(dble(na)*mass)
      nax   = maxval(na)
      allocate(tau(3,nax,nsp))
      allocate(taui(3,nax,nsp))
      allocate(sqd(nax,nsp))
      allocate(msqd(nsp))

      open(unit=34,file=atofile,status='old')

      tau = 0.0d0
      taui = 0.0d0

      do t = 1, nskip
        read(34,'(a)') dummy
        do is = 1, nsp
          do ia = 1, na(is)
            read(34,*) (tau(k,ia,is),k=1,3) 
          end do
        end do
      end do

      do t = nskip+1, nstep

        read(34,'(a)') dummy
        do is = 1, nsp
          do ia = 1, na(is)
            read(34,*) (tau(k,ia,is),k=1,3) 
          end do
        end do

        cdm = 0.d0
        do is = 1, nsp
          do ia = 1, na(is)
            cdm(:) = cdm(:) + tau(:,ia,is)*mass(is)
            !cdm(:) = cdm(:) + tau(:,ia,is)
          end do
        end do
        cdm = cdm / masst

        if(t.eq.(nskip+1)) then
          cdmi = cdm
          do is = 1, nsp
            do ia = 1, na(is)
              taui(:,ia,is) = tau(:,ia,is) - cdmi(:)
            end do
          end do
        end if

        time = dble(t) * dt + toffset
        sqd  = 0.D0
        msqd = 0.D0
        do is = 1, nsp
          do ia = 1, na(is)
            r(:) = tau(:,ia,is) - cdm(:)
            sqd(ia,is)  = sum( ((r(:)-taui(:,ia,is)))**2 )
          end do
          msqd(is) = sum(sqd(:,is))/(dble(na(is)))
        end do
        write(4,*) time,time*picosecond,(msqd(is)*(0.52917)**2,is=1,nsp)
      end do
      close(unit=34)
      deallocate(na)
      deallocate(mass)
      deallocate(tau)
      deallocate(taui)
      deallocate(sqd)
      deallocate(msqd)

 100  format(f9.1,2x,f9.7,8f12.6)

      END program diff

