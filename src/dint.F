      program dint
c     DiNT:  Adiabatix and Nonadiabatic Trajectories
c     Main program.  This routine calls HEADER and READIN once, then
c     calls PREMOL once per AG to precalculate some constants.
c     DRIVER is called NTRAJ times, once for each trajectory.
c     FINALSTATE is called after each trajectory to write some data.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
      include 'c_traj.f'
#include <sprng_f.h>

      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      integer im,initrng,minitrng,i,ia
      double precision t1,t2,t3,cputime,t3old,r

ccccc MPI
      call MPI_INIT(ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

      IF (my_id.eq.0) then
        write(6,"(a, i3)") " MPI num procs: ", nproc
!        write(6,"(a, i3)") " MPI current proc: ", my_id
        call header
        write(6,'("...entering READIN...",/)')
      ENDIF

      call readin
      IF (my_id.eq.0) then
        write(6,'("...exiting READIN...",/)')
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      IF (my_id.eq.0) then
        call timing(t1)
        write(6,'("...entering PREMOL...",/)')
      ENDIF

c     prepare each molecule
c     compute quantities that are the same for all trajectories
c     initialize electronic coordinates
      call initelec
!      print *,"nmol = ",nmol," according to proc",my_id
      do im = 1, nmol
        call premol(im)
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      IF (my_id.eq.0) then
        write(6,'("...exiting PREMOL...",/)')
        call timing(t2)
        write(6,'(" CPU time in PREMOL is ",f10.5," s",/)')(t2-t1)
      ENDIF

c     redirect for N(E)
      if (methflag.eq.10) then
        initrng=1
        rng_stream = init_sprng(initrng,maxtraj,ranseed,SPRNG_DEFAULT)
        im=1
        call mc(im)
        go to 999
      endif 
      
      t3 = 0.d0
      cputime = 0.d0
c     loop over trajectories
      nistep=0
      nistepw=0

      DO itraj=1,ntraj
        IF (my_id.eq.0) then
        write(6,*)"***TRAJECTORY ",itraj," OUT OF ",ntraj,"***"
        write(6,*)
        ENDIF
        if (tflag(2).eq.1) then
          IF (my_id.eq.0) then
          write(6,*)"Restarting trajectory #",trajlist(itraj)
          write(6,*)
          ENDIF
          initrng = trajlist(itraj)-1
        else
          initrng = itraj-1
        endif

c       initialize this trajectory
        nsurf = nsurf0
        hstep = hstep0
c       initialize electronic coordinates
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call initelec
c       initialize the randum number generator
        rng_stream = init_sprng(initrng,maxtraj,ranseed,SPRNG_DEFAULT)

        IF (my_id.eq.0) write(6,'("...entering INITMOL...",/)')
        do im = 1 , nmol
c         generate specific initial coordinates for each molecule
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call initmol(im)
        enddo
        IF (my_id.eq.0) then
        write(6,'("...exiting INITMOL...",/)')

        write(6,'("...entering ORIENT...",/)')
        ENDIF
        do im = 1 , nmol
c         orient AGs with respect to each other
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call orient(im)
        enddo
        IF (my_id.eq.0) then
        write(6,'("...exiting ORIENT...",/)')

        write(6,'("...entering PREPTRAJ...",/)')
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call preptraj
        IF (my_id.eq.0) then
        write(6,'("...exiting PREPTRAJ...",/)')

c       propagate this trajectory
        write(6,'("...entering DRIVER...",/)')
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call driver
        IF (my_id.eq.0) then
        write(6,'("...exiting DRIVER...",/)')
        t3old = t3
        call timing(t3)
        cputime = cputime + (t3-t3old)
        write(6,'(" CPU time for this trajectory is ",f10.5," s",/)')
     &    (t3-t3old)
        ENDIF

c       analyze trajectory
        IF (my_id.eq.0) write(6,'("...entering FINALSTATE...",/)')
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call finalstate
        IF (my_id.eq.0) write(6,'("...exiting FINALSTATE...",/)')

c     end loop over trajectories
      enddo

c     from ransamp.F

c MPI
      IF (my_id.eq.0) THEN
      write(6,'(" CPU time for all trajectories is ",f10.5," s or ",
     & f10.5," s per trajectory.",/)')
     & cputime,cputime/dble(ntraj)
      ENDIF

 999  continue

c MPI
      IF (my_id.eq.0) write(6,'("Well done!")')

ccccc MPI
      call MPI_FINALIZE(ierr)

      end
