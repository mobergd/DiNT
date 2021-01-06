      subroutine preptraj

c Put center of mass at origin

      implicit none
      include 'param.f'
      include 'c_traj.f'
      include 'c_sys.f'
 
c MPI
      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      double precision xtot,ytot,ztot,mtot

      integer i

ccccc MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

      call MPI_BCAST(pp, 3*nat, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call gettemp(pp,mm,nat,temp,ke)
      IF (my_id.eq.0) THEN
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF

      IF (my_id.eq.0) THEN
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,nat
        xtot = xtot + mm(i)*xx(1,i)
        ytot = ytot + mm(i)*xx(2,i)
        ztot = ztot + mm(i)*xx(3,i)
        mtot = mtot + mm(i)
      enddo
      do i=1,nat
        xx(1,i) = xx(1,i) - xtot/mtot
        xx(2,i) = xx(2,i) - ytot/mtot
        xx(3,i) = xx(3,i) - ztot/mtot
      enddo
      ENDIF
      call MPI_BCAST(xx, 3*nat, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call gettemp(pp,mm,nat,temp,ke)
      IF (my_id.eq.0) THEN
      write(6,*)"Placed CoM at origin"
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF

c     remove overall momenta
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,nat
        xtot=xtot+pp(1,i)
        ytot=ytot+pp(2,i)
        ztot=ztot+pp(3,i)
        mtot = mtot + mm(i)
      enddo
      do i=1,nat
        pp(1,i)=pp(1,i)-xtot*mm(i)/mtot
        pp(2,i)=pp(2,i)-ytot*mm(i)/mtot
        pp(3,i)=pp(3,i)-ztot*mm(i)/mtot
      enddo
      call gettemp(pp,mm,nat,temp,ke)
      IF (my_id.eq.0) THEN
      write(6,*)"Removed CoM motion"
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF

      return

      end
