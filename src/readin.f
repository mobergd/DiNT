      subroutine readin

c This subroutine reads from the standard input.  This is the only 
c subroutine that reads data.  Data are stored for later use in
c common blocks.

      implicit none

      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
      include 'c_traj.f'
      include 'c_output.f'

c MPI
      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

c local variables
      integer i,ii,im,k,ntmp,itmp(99),nrecl
      double precision rdum,murel,prel,dum
      character*2 chardum

      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

c read potential interface flag
      IF (my_id.eq.0) THEN
      write(6,*)"Potential interface"
      write(6,*)"-------------------"
      read(5,*)potflag
      write(6,*)"POTFLAG = ",potflag
      ENDIF
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(potflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      if (potflag.eq.0) then
        IF (my_id.eq.0) write(6,*)"Using the POTLIB MM-HO-1 interface"
      elseif (potflag.eq.2) then
        IF (my_id.eq.0) write(6,*)"Using the POTLIB MM-HE-1 interface"
      elseif (potflag.eq.1) then
        IF (my_id.eq.0) write(6,*)"Using the POTLIB 3-2V interface"
      elseif (potflag.eq.3) then
        IF (my_id.eq.0) write(6,*)
     &    "Using the multistate diabatic MM-HE-1 interface"
      elseif (potflag.eq.4) then
        IF (my_id.eq.0) write(6,*)
     &    "Using the multistate adiabatic MM-HE-1 interface"
c     ADD YOUR OWN INTERFACE HERE, USING POTFLAG = -1
c     YOU NEED ONLY MODIFY GETPEM
c      elseif (potflag.eq.-1) then
c      write(6,*)"Using a user-designed interface"
      else
        IF (my_id.eq.0) write(6,*)
     &    "POTFLAG = ",potflag," is not supported"
        stop
      endif
      IF (my_id.eq.0) write(6,*)

!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call prepot

c LDOFRAG
      ldofrag = .false.
      ldofrag = .true.

c initial surface, number of coupled surfaces, and electronic representation flag
      IF (my_id.eq.0) THEN
      read(5,*)nsurf0,nsurft,methflag,repflag
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nsurf0, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nsurft, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(methflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(repflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      IF (my_id.eq.0) THEN
      write(6,*)"Electronic surfaces"
      write(6,*)"-------------------"
      write(6,*)"Initial electronic surface = ",nsurf0
      write(6,*)"Number of electronic surfaces = ",nsurft
      write(6,*)"REPFLAG = ",repflag
      ENDIF

      if (repflag.eq.0) then
        IF (my_id.eq.0)
     &    write(6,*)"The adiabatic representation will be used."
      else if (repflag.eq.1) then
        IF (my_id.eq.0)
     &    write(6,*)"The diabatic representation will be used."
        if (potflag.eq.4) then
          IF (my_id.eq.0)
     &      write(6,*)"The diabatic representation cannot be",
     &      " used with POTFLAG=4"
          stop
        endif
      else
        IF (my_id.eq.0) 
     &    write(6,*)"REPFLAG = ",repflag," is not supported"
        stop
      endif

      if (methflag.eq.0) then
        IF (my_id.eq.0)
     &    write(6,*)"METHFLAG = 0:  Single surface propagation"
      elseif (methflag.eq.1) then
        IF (my_id.eq.0) THEN
          write(6,*)"METHFLAG = 1:  Electronically nonadiabatic",
     &    " calculation"
          write(6,*)"Tully's fewest switches method"
        ENDIF
      elseif (methflag.eq.2) then
        IF (my_id.eq.0) THEN
          write(6,*)"METHFLAG = 2:  Electronically nonadiabatic",
     &    " calculation"
          write(6,*)"Semiclassical Ehrenfest method"
        ENDIF
      elseif (methflag.eq.3) then
        IF (my_id.eq.0) THEN
          write(6,*)"METHFLAG = 3:  Electronically nonadiabatic",
     &    " calculation"
          write(6,*)"Self consistent decay of mixing method"
        ENDIF
      elseif (methflag.eq.4) then
        IF (my_id.eq.0) THEN
          write(6,*)"METHFLAG = 4:  Electronically nonadiabatic",
     &    " calculation"
          write(6,*)"Coherent switching with decay of mixing method"
        ENDIF
      elseif (methflag.eq.5) then
        IF (my_id.eq.0) THEN
        write(6,*)"METHFLAG = 5:  Electronically nonadiabatic",
     &  " calculation"
        write(6,*)"Fewest switches with time uncertainty method"
        ENDIF
      elseif (methflag.eq.10) then
        IF (my_id.eq.0)
     &    write(6,*)"METHFLAG = 10: Monte Carlo sampling of W(E)"
      else
        IF (my_id.eq.0)
     &    write(6,*)"METHFLAG = ",methflag," is not allowed"
        stop
      endif
      if (nsurf0.gt.nsurft) then
        IF (my_id.eq.0)
     &    write(6,*)"Initial surface greater than total number!"
        stop
      endif
      if (nsurft.gt.mnsurf) then
        IF (my_id.eq.0)
     &    write(6,*)"Number of surfaces > MNSURF in PARAMS"
        stop
      endif

      IF (my_id.eq.0) THEN
      if (methflag.eq.10) then
        read (5,*)mccurv
        if (mccurv.gt.0) then
          do i=1,mccurv
            read(5,*)mcmode(i),mctype(i),mcpar(1,mcmode(i)),
     &               mcpar(2,mcmode(i))
            if (mctype(i).eq.1) then
             write(6,*)"Treating mode ",mcmode(i)," as a torsion about",
     & " atoms",mcpar(1,mcmode(i))," and ",mcpar(2,mcmode(i))
            elseif (mctype(i).eq.-1) then
              write(6,*)"Skipping mode ",mcmode(i),"."
            endif
          enddo
        endif
      endif
      write(6,*)
      ENDIF

!      do i=1,nsurft*2  ! why is this here? This should be in INITMOL
!        cre(i) = 0.d0
!        cim(i) = 0.d0
!      enddo
!      cre(nsurf0) = 1.d0

c read integration parameters
      IF (my_id.eq.0) THEN
      write(6,*)"Integrator information"
      write(6,*)"----------------------"
      read(5,*)intflag
      write(6,*)"INTFLAG = ",intflag
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(intflag, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

      if (intflag.eq.1) then
        IF (my_id.eq.0) THEN
        read(5,*)hstep0,nprint
        write(6,*)"Using the 4th-order Runge-Kutta integrator",
     &  " with "
        write(6,*)"HSTEP=",hstep0," fs"
c       convert timestep to au
        hstep0 = hstep0/autofs
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(hstep0, 1, MPI_DOUBLE_PRECISION,
     &                         0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nprint, 1, MPI_INTEGER,
     &                         0, MPI_COMM_WORLD, ierr)
      elseif (intflag.eq.0) then
        IF (my_id.eq.0) THEN
        read(5,*)hstep0,eps,nprint
        write(6,*)"Using the Bulirsch-Stoer integrator",
     &  " with "
        write(6,*)"initial HSTEP=",hstep0," fs"
        write(6,*)"and tolerence EPS=",eps
c       convert timestep to au
        hstep0 = hstep0/autofs
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!        IF (my_id.ne.0) THEN
        call MPI_BCAST(hstep0, 1, MPI_DOUBLE_PRECISION,
     &                         0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(eps, 1, MPI_DOUBLE_PRECISION,
     &                      0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nprint, 1, MPI_INTEGER,
     &                         0, MPI_COMM_WORLD, ierr)
!        ENDIF
      else
        IF (my_id.eq.0)
     &    write(6,*)"INTFLAG = ",intflag," is not supported"
        stop
      endif
      IF (my_id.eq.0) THEN  
      write(6,*)"Write information every NPRINT=",nprint," steps"
      write(6,*)
      ENDIF

c read random number seed
      IF (my_id.eq.0) THEN  
      write(6,*)"RNG information"
      write(6,*)"---------------"
      read(5,*)ranseed
      write(6,*)"Random number seed is ",ranseed
      write(6,*)
      ENDIF

      call MPI_BCAST(ranseed, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD, ierr)

c read trajectory conrol flag
      IF (my_id.eq.0) THEN  
      write(6,*)"Trajectory control flags"
      write(6,*)"-----------------------"
      read(5,*)ntraj,(tflag(i),i=1,4)
      write(6,*)"Running ",ntraj," trajectories"
      ENDIF
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ntraj, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(tflag, 4, MPI_INTEGER,
     &               0, MPI_COMM_WORLD, ierr)
      if (ntraj.gt.mntraj) then
        IF (my_id.eq.0) write(6,*)"NTRAJ = ",ntraj," > MNTRAJ = ",mntraj
        stop
      endif

      if (tflag(1).eq.0) then
      elseif (tflag(1).eq.1) then
        IF (my_id.eq.0)
     &  write(6,*)"TFLAG(1) = 1:  Momenta set to zero at every step,",
     &  " i.e., steepest descent minimization"
      elseif (tflag(1).eq.2) then
        IF (my_id.eq.0) write(6,*)"TFLAG(1) = 2:  Temperature ramping"
      elseif (tflag(1).eq.3) then
        IF (my_id.eq.0) write(6,*)"TFLAG(1) = 3:  Andersen thermostat"
c      elseif (tflag(1).eq.4) then
c        write(6,*)"TFLAG(1) = 4:  MSX finder"   ! Doesn't currently work
      else
        IF (my_id.eq.0) write(6,*)"TFLAG(1) != 0-3 is not supported"
        stop
      endif

      if (tflag(2).eq.0) then
      elseif (tflag(2).eq.1) then
        IF (my_id.eq.0) write(6,*)"TFLAG(2) = 1:  Restart trajectories"
      else
        IF (my_id.eq.0) write(6,*)"TFLAG(2) != 0 or 1 is not supported"
        stop
      endif
      if (tflag(3).eq.0) then
      elseif (tflag(3).eq.1) then
        IF (my_id.eq.0)
     &    write(6,*)"TFLAG(3) = 1:  Photoexcited trajectory"
      else
        IF (my_id.eq.0) write(6,*)"TFLAG(3) != 0 or 1 is not supported"
        stop
      endif

      if (tflag(1).eq.2) then
        IF (my_id.eq.0) THEN
        read(5,*)ramptime,rampfact,nramp
        write(6,*)"The momentum will be ramped every RAMPTIME = ",
     &   ramptime," fs by "
        write(6,*)"a factor of RAMPFACT = ",rampfact,
     &  " a total of NRAMP = ",nramp," times."
        ramptime=ramptime/autofs
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ramptime, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(rampfact, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(nramp, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD, ierr)      
      endif

      if (tflag(1).eq.3) then
        IF (my_id.eq.0) THEN
        read(5,*)andersen_temp,andersen_freq,scandth
        write(6,*)"Temp will be thermostated to temp = ",andersen_temp
        write(6,*)"Andersen frequency parameter = ",
     &   andersen_freq," 1/fs"
c       convert to 1/a.u.time
        andersen_freq=andersen_freq*autofs
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(andersen_temp, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(andersen_freq, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(scandth, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
        if (scandth.gt.0.d0) then
          IF (my_id.eq.0)
     &      write(6,101)"Resulting momenta scaled such that the total",
     &      " energy = ",scandth," eV plus the initial potential energy"
          scandth=scandth/autoev
        else
          IF (my_id.eq.0) write(6,*)"Resulting momenta are not scaled"
        endif
      endif
      if (tflag(1).eq.4) then
        IF (my_id.eq.0) THEN
        write(6,*)"Gradient will be modified to point toward MSX"
        write(6,*)"Momentum will be zeroed at every step"
        ENDIF
        if (nsurft.ne.2) then
          IF (my_id.eq.0) write(6,*)"NSURFT must be 2 for TFLAG(1) = 4"
          stop
        endif
      endif

      if (tflag(2).eq.1) then
        IF (my_id.eq.0) THEN
        read(5,*)(trajlist(k),k=1,ntraj)
        write(6,*)"Restarting trajectories: ",(trajlist(k),k=1,ntraj)
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(trajlist, ntraj, MPI_INTEGER,
     &               0, MPI_COMM_WORLD, ierr)      
        maxtraj = trajlist(1)
        do k=1,ntraj
          maxtraj = max0(maxtraj,trajlist(k))
        enddo
      else
        maxtraj = ntraj
      endif
      IF (my_id.eq.0) write(6,*)

      if (tflag(3).eq.1) then
        IF (my_id.eq.0) THEN
        read(5,*)ntarget,ephoton,wphoton
        write(6,*)"Exciting every trajectory from initial state ",
     &   nsurf0," to state ",ntarget,
     &  " with a photon energy of ",ephoton," +/- ",wphoton," eV"
c       convert eV to hartree
        ephoton=ephoton/autoev
        wphoton=wphoton/autoev
        ENDIF
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ntarget, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(ephoton, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(wphoton, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
!        ENDIF
      endif

      stodecoflag=.false.
      if (tflag(4).eq.1) then
        IF (my_id.eq.0) write(6,*)
     &  "Using the stochastic decoherence (SD) method"
        stodecoflag=.true.
      elseif (tflag(4).lt.0) then
        IF (my_id.eq.0) write(6,*)
     &  "Gathering seam crossing statistics for the first ",
     & -tflag(4)," seam crossings."
        if(repflag.ne.1) then
          IF (my_id.eq.0) write(6,*)
     &   "The diabatic representation must be used for TFLAG(4) < 1"
          stop
        endif
      endif

      IF (my_id.eq.0) THEN
      write(6,*)
      write(6,*)"Atom group (AG) specification"
      write(6,*)"-----------------------------"
c read number of molecules (atom groups)
      read(5,*)nmol,ezero
      write(6,*)"This is a ",nmol," AG calculation"
      write(6,*)"The zero of energy for the supermolecule is ",
     & ezero," eV"
      ezero=ezero/autoev
      write(6,*)
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nmol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ezero, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)

      if (nmol.gt.mnmol) then
       IF (my_id.eq.0) write(6,*)"ERROR:  nmol > mnmol, ",nmol,">",mnmol
       stop
      endif
      if (nmol.le.0) then
        IF (my_id.eq.0) write(6,*)"ERROR:  nmol must be > 0, nmol=",nmol
        stop
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
c read info for each molecule
      nat = 0
      DO im=1,nmol
        IF (my_id.eq.0) THEN
        write(6,*)"AG ",im
        write(6,*)"-------"
        read(5,*)natom(im),initx(im),initp(im),initj(im),ezeroim(im)
        ENDIF
!        natom1=natom(im)
!        initx1=initx(im)
!        initp1=initp(im)
!        initj1=initj(im)
!        ezero1=ezeroim(im)

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(natom, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(initx, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(initp, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(initj, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(ezeroim, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      

        iatom(im)=nat
        nat=nat+natom(im)
        IF (my_id.eq.0) write(6,*)
     &  "This AG has ",natom(im)," atoms (",
     &  "atoms ",iatom(im)+1," to ",nat,")"
        if (natom(im).le.0) then
          IF (my_id.eq.0) write(6,*)
     &    "ERROR: cant have natom < 1, natom(",im,")=",natom(im)
          stop
        endif
        if (nat.gt.mnat) then
          IF (my_id.eq.0) write(6,*)"ERROR:  nat > mnat, ",nat,">",mnat
          stop
        endif
        IF (my_id.eq.0) write(6,*)
     &  "Zero of energy for this AG is ",ezeroim(im)," eV"
        ezeroim(im)=ezeroim(im)/autoev

        IF (my_id.eq.0) write(6,*)
     &   "Initial conditions prepared using INITx = ",
     &   initx(im)," INITp = ",initp(im)," and INITj = ",initj(im)

c INITx controls the initial coordinates
c read different info depending on the value of INITx(IM)
        if (initx(im).eq.0) then
c         crystal structure
          IF (my_id.eq.0) THEN
          write(6,*)"INITx = 0:  Reading in initial coordinates:"
          write(6,102)"    # symb","mass (amu)","x (A)","y (A)","z (A)"
          ENDIF
 102      format(a10,4a12)
          mmag(im)=0.d0
          IF (my_id.eq.0) THEN
          do i=1,natom(im)
            ii=iatom(im)+i
c           read in atomic symbol, coords, and mass
            read(5,*)symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
            write(6,103)ii,symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
 103        format(i5,a5,4f12.5)
c           convert from amu to au
            mm(ii)=mm(ii)*amutoau
c           compute total mass for this AG
            mmag(im)=mmag(im)+mm(ii)
c           convert from A to bohr
            do k=1,3
              xx0(k,ii)=xx0(k,ii)/autoang
            enddo
          enddo
          write(6,*)
          ENDIF
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
!          call MPI_BCAST(ii, 1, MPI_INTEGER,
!     &                 0, MPI_COMM_WORLD, ierr)      
!          do i=1,natom(im)
!          call MPI_BCAST(symbol, natom(im)*2, MPI_CHARACTER,
!     &                   0, MPI_COMM_WORLD, ierr)
!          enddo
          call MPI_BCAST(mm, natom(im), MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(mmag, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(xx0, 3*natom(im), MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
        elseif (initx(im).eq.1) then
c       random spherical clusters
          IF (my_id.eq.0) read(5,*) rdum
          call MPI_BCAST(rdum, 1, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          rran(1) = rdum*dble(natom(im))**(1.d0/3.d0)
          rran(2) = rdum*0.7d0
          rran(3) = rdum*1.3d0
          IF (my_id.eq.0) THEN
          write(6,*)"INITx = 1:  Generate random spherical clusters",
     &  " of radius ",rran(1)," A"
          write(6,*)"Include ",rran(2)," < rmin < ",rran(3)," A"
          write(6,102)"    # symb","mass (amu)"
          ENDIF
c         convert from A to bohr
          do i=1,3
            rran(i)=rran(i)/autoang
          enddo
          mmag(im)=0.d0
          IF (my_id.eq.0) THEN
          do i=1,natom(im)
            ii=iatom(im)+i
c           read in atomic symbol, and mass
            read(5,*)symbol(ii),mm(ii)
            write(6,103)ii,symbol(ii),mm(ii)
c           convert from amu to au
            mm(ii)=mm(ii)*amutoau
c           compute total mass for this AG
            mmag(im)=mmag(im)+mm(ii)
          enddo
          write(6,*)
          ENDIF
          call MPI_BCAST(ii, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(symbol, ii, MPI_CHARACTER,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(mm, ii, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(mmag, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        elseif (initx(im).eq.2) then
c         normal modes
          IF (my_id.eq.0) write(6,*)
     &    "INITx = 2: State selected initial conditions"
          ntmp = 3*natom(im)-6
          if (natom(im).eq.2) ntmp = 1
c          read(5,*)lreadhess,nmtype(im),(nmqn(k,im),k=1,ntmp),scale0im(im)
          IF (my_id.eq.0)
     &    read(5,*)lreadhess,nmtype(im),(nmqn(k,im),k=1,ntmp)
          call MPI_BCAST(lreadhess, 1, MPI_LOGICAL,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(nmtype, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(nmqn, ntmp*nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
c          if (scale0im(im).gt.0.d0) then
c            write(6,101)"Resulting momenta scaled such that the total",
c       &     " energy = ",scale0im(im)," eV"
c            scale0im(im)=scale0im(im)/autoev
c          else
c            write(6,*)"Resulting momenta are not scaled"
c          endif
          IF (my_id.eq.0) THEN
          if (nmtype(im).eq.1) then
            write(6,*)"Initial geometry is a saddle point"
            write(6,*)"Unbound mode will be assigned an energy of ",
     &      nmqn(ntmp,im)," eV"
            nmqn(ntmp,im)=nmqn(ntmp,im)/autoev
            write(6,104)"Initial quantum numbers for the ",ntmp-1,
     &      " bound modes are:",(nmqn(k,im),k=1,ntmp-1)
          else
            write(6,*)"Initial geometry is a local minimum"
            write(6,104)"Initial quantum numbers for the ",ntmp,
     &      " bound modes are:",(nmqn(k,im),k=1,ntmp)
          endif
          write(6,*)
          ENDIF
 104      format(1x,a,i5,a,1000f8.2)
          mmag(im)=0.d0
          IF (my_id.eq.0) THEN
          do i=1,natom(im)
            ii=iatom(im)+i
c           read in atomic symbol, coords, and mass
            read(5,*)symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
            write(6,113)ii,symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
 113        format(i5,a5,4f12.5)
c           convert from amu to au
            mm(ii)=mm(ii)*amutoau
c           compute total mass for this AG
            mmag(im)=mmag(im)+mm(ii)
c           convert from A to bohr
            do k=1,3
              xx0(k,ii)=xx0(k,ii)/autoang
            enddo
          enddo
          write(6,*)"For INITx = 2, INITp is not used and is set to -1"
          write(6,*)
          ENDIF
          call MPI_BCAST(ii, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(symbol, ii, MPI_CHARACTER,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(mm, ii, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(mmag, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(xx0, 3*ii, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          initp(im) = -1
      elseif (initx(im).eq.3) then
c       atom-diatom initial conditions
        IF (my_id.eq.0) 
     &   write(6,*)"INITx = 3:  Atom-diatom initial conditions"
        if (natom(im).ne.3) then
          IF (my_id.eq.0) write(6,*)
     &     "INITx = 3 must be used with a triatomic molecule"
          stop
        endif
        if (nmol.ne.1) then
          IF (my_id.eq.0) write(6,*)
     &     "INITx = 3 may not be used with multiple AGs"
          stop
        endif
        mmag(im)=0.d0
        IF (my_id.eq.0) THEN
        write(6,102)"    # symb","mass (amu)"
        do i=1,natom(im)
          ii=iatom(im)+i
c         read in atomic symbol, and mass
          read(5,*)symbol(ii),mm(ii)
          write(6,103)ii,symbol(ii),mm(ii)
c         convert from amu to au
          mm(ii)=mm(ii)*amutoau
c         compute total mass for this AG
          mmag(im)=mmag(im)+mm(ii)
        enddo
        read(5,*)escatad,vvad,jjad,rrad,arrad
        write(6,*)"Total energy = ",escatad," eV"
        write(6,*)"Initial v = ",vvad
        write(6,*)"Initial j = ",jjad
        write(6,*)"Initial distance = ",rrad
        write(6,*)"Initial arrangement = ",arrad
        write(6,*)"For INITx = 3, INITp is not used and is set to -1"
        escatad=escatad/autoev
        rrad=rrad/autoang
        write(6,*)
        ENDIF
        call MPI_BCAST(ii, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(symbol, ii, MPI_CHARACTER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(mm, ii, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(mmag, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(escatad, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(vvad, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(jjad, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(rrad, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(arrad, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        initp(im) = -1
      elseif (initx(im).eq.4) then
c       NOT IMPLEMENTED. DO NOT USE.
c       diatom initial conditions
        write(6,*)"INITx = 4:  Diatom initial conditions"
        if (natom(im).ne.2) then
          write(6,*)"INITx = 4 must be used with a diatomic molecule"
          stop
        endif
        read(5,*)vvdi(im),jjdi(im),rmindi(im)
        if (vvdi(im).ge.0.d0) then
          write(6,*)"Initial v = ",vvdi(im)
        else
          write(6,*)"Initial v selected from thermal distribution, T = "
     &        ,-vvdi(im)," K"
        endif
        if (jjdi(im).ge.0.d0) then
          write(6,*)"Initial j = ",jjdi(im)
        else
          write(6,*)"Initial j selected from thermal distribution, T = "
     &        ,-jjdi(im)," K"
        endif
        write(6,*)"Estimated R = ",rmindi(im)
        rmindi(im)=rmindi(im)/autoang
        write(6,*)
        write(6,102)"    # symb","mass (amu)"
        mmag(im)=0.d0
        do i=1,natom(im)
          ii=iatom(im)+i
c         read in atomic symbol, and mass
          read(5,*)symbol(ii),mm(ii)
          write(6,103)ii,symbol(ii),mm(ii)
c         convert from amu to au
          mm(ii)=mm(ii)*amutoau
c         compute total mass for this AG
          mmag(im)=mmag(im)+mm(ii)
        enddo
      elseif (initx(im).eq.5) then
c       normal modes
        IF (my_id.eq.0) write(6,*)
     &   "INITx = 5: Vibrational thermal state selected ",
     &   "initial conditions"
        ntmp = 3*natom(im)-6
        do k=1,ntmp
          nmqn(k,im)=0.d0
        enddo
c        read(5,*)lreadhess,temp0im(im),scale0im(im),sampwell(im),
        sampwell(im)=1
        IF (my_id.eq.0) THEN
        read(5,*) lreadhess,nmtype(im),temp0im(im),scale0im(im),lems(im)
        if (temp0im(im).lt.0.d0) then
          write(6,101)"Resulting momenta scaled such that the ",
     &     " harmonic energy = ",scale0im(im)," eV relative to ",
     &     " the bottom of the well"
          scale0im(im)=scale0im(im)/autoev
        else
          write(6,107)"Target temperature is ",temp0im(im)," K"
        endif
        if (lems(im)) then
          read(5,*)emicr(im),ninc(im),nbrea(im),nemstot(im)
          write(6,*)"Efficient Microcanonical Sampling (EMS) ",
     &     "will be used to refine the resulting coordinates ",
     &     "and momenta to generate a microcanonical ensemble ",
     &     "with Etot = ",emicr(im)," eV"
          emicr(im)=emicr(im)/autoev
          emsw(im)=0.d0
        endif
        write(6,*)
        write(6,*)"Geometry of the stationary point"
        write(6,102)"    # symb","mass (amu)","x (A)","y (A)","z (A)"
        ENDIF
 107    format(1x,a,f15.5,a)
        call MPI_BCAST(lreadhess, 1, MPI_LOGICAL,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(nmtype, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(temp0im, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(scale0im, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(emicr, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(emsw, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(nbrea, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(nemstot, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(lems, nmol, MPI_LOGICAL,
     &                 0, MPI_COMM_WORLD, ierr)      
        mmag(im)=0.d0
        IF (my_id.eq.0) THEN
        do i=1,natom(im)
          ii=iatom(im)+i
c         read in atomic symbol, coords, and mass
          read(5,*)symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
          write(6,113)ii,symbol(ii),mm(ii),(xx0(k,ii),k=1,3)
c         convert from amu to au           
          mm(ii)=mm(ii)*amutoau
c         compute total mass for this AG
          mmag(im)=mmag(im)+mm(ii)
c         convert from A to bohr
          do k=1,3
            xx0(k,ii)=xx0(k,ii)/autoang
          enddo
c          write(6,*)
        enddo
        ENDIF
        call MPI_BCAST(ii, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(symbol, ii, MPI_CHARACTER,
     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(mm, ii, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(mmag, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(xx0, 3*ii, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)
        if (sampwell(im).eq.2) then
          IF (my_id.eq.0) read(5,*)relwell2(im),emin2(im)
          call MPI_BCAST(relwell2, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(emin2, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          emin2(im)=emin2(im)/autoev
          relwell1(im)=1.d0-relwell2(im)
          IF (my_id.eq.0) THEN
          write(6,*)"Second well"
          do i=1,natom(im)
            ii=iatom(im)+i
            read(5,*)(xx02(k,ii),k=1,3)
            write(6,113)ii,symbol(ii),mm(ii)/amutoau,(xx02(k,ii),k=1,3)
            do k=1,3
              xx02(k,ii)=xx02(k,ii)/autoang
            enddo
          enddo
          ENDIF
          call MPI_BCAST(ii, 1, MPI_INTEGER,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(xx02, 3*ii, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)
        endif
        IF (my_id.eq.0) THEN
        write(6,*)
        write(6,*)"For INITx = 5, INITp is not used and is set to -1"
        write(6,*)
        ENDIF
        initp(im) = -1
      elseif (initx(im).eq.6) then
        IF (my_id.eq.0) THEN
        write(6,*)"INITx = 6: Read in initial coordinates and momenta ",
     &   "from a separate files"
        read(5,*)samptot(im),lbinsamp(im),sampfilexx(im),sampfilepp(im),
     &   lems(im)
        write(6,*)"Sampling from ",samptot(im)," data contained in ",
     &    "the XYZ file ",sampfilexx(im)," and the PxPyPz file ",
     &    sampfilepp(im)
        ENDIF
c     &   sampwell(im),lems(im)
        call MPI_BCAST(samptot, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(lbinsamp, nmol, MPI_LOGICAL,
     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampfilexx, nmol, MPI_CHARACTER,
     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampfilepp, nmol, MPI_CHARACTER,
     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(lems, nmol, MPI_LOGICAL,
     &                   0, MPI_COMM_WORLD, ierr)      
        sampwell(im)=1
        if (sampwell(im).eq.2) then 
          IF (my_id.eq.0) 
     &     read(5,*)samptot2(im),lbinsamp2(im),sampfilexx2(im),
     &     sampfilepp2(im),relwell2(im)
          call MPI_BCAST(samptot2, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(lbinsamp2, nmol, MPI_LOGICAL,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(sampfilexx2, nmol, MPI_CHARACTER,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(sampfilepp2, nmol, MPI_CHARACTER,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(relwell2, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          relwell1(im)=1.d0-relwell2(im)
        endif
        if (lems(im)) then
          IF (my_id.eq.0) THEN
          read(5,*)emicr(im),ninc(im),nbrea(im),nemstot(im)
          write(6,*)"Efficient Microcanonical Sampling (EMS) ",
     &     "will be used to refine the resulting coordinates ",
     &     "and momenta to generate a microcanonical ensemble ",
     &     "with Etot = ",emicr(im)," eV"
          emicr(im)=emicr(im)/autoev
          ENDIF
          emsw(im)=0.d0
          call MPI_BCAST(emicr, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(ninc, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(nbrea, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(nemstot, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        endif
        IF (my_id.eq.0) THEN
        write(6,*)
        write(6,102)"    # symb","mass (amu)"
        mmag(im)=0.d0
        do i=1,natom(im)
          ii=iatom(im)+i
c         read in atomic symbol, coords, and mass
          read(5,*)symbol(ii),mm(ii)
          write(6,103)ii,symbol(ii),mm(ii)
c         convert from amu to au
          mm(ii)=mm(ii)*amutoau
c         compute total mass for this AG
          mmag(im)=mmag(im)+mm(ii)
c         convert from A to bohr
        enddo
        write(6,*)
        write(6,*)"For INITx = 6, INITp is not used and is set to -1"
        write(6,*)
        ENDIF
!        call MPI_BCAST(ii, 1, MPI_INTEGER,
!     &                 0, MPI_COMM_WORLD, ierr)      
!        call MPI_BCAST(symbol, ii, MPI_CHARACTER,
!     &                   0, MPI_COMM_WORLD, ierr)      
!        call MPI_BCAST(mm, ii, MPI_DOUBLE_PRECISION,
!     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(mmag, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
        initp(im) = -1
      else
        IF (my_id.eq.0) write(6,*)
     &     "ERROR:  initial conditions flag INITx = ",initx(im),
     &      " is not supported"
        stop
      endif

c INITp controls the initial momenta
      if (initp(im).eq.0) then
c       initial target temperature
        IF (my_id.eq.0) THEN
        read(5,*)temp0im(im),escale0im(im)
        write(6,101)"INITp = 0:  Target temp for random thermal ",
     &    "dist is ",temp0im(im),"K"
        if (escale0im(im).gt.0.d0) then
          write(6,101)"Resulting momenta scaled such that the ",
     &     "kinetic energy = ",escale0im(im)," eV"
          escale0im(im)=escale0im(im)/autoev
        else
          write(6,*)"Resulting momenta are not scaled"
        endif
        write(6,*)
        ENDIF
        call MPI_BCAST(temp0im, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(escale0im, nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
 101    format(1x,a,a,f12.5,1x,a)
      elseif (initp(im).eq.1) then
        IF (my_id.eq.0) THEN
        write(6,*)"INITp = 1:  Initial momenta set to zero"
        write(6,*)
        ENDIF
      elseif (initp(im).eq.2) then
        IF (my_id.eq.0) THEN
        write(6,*)"INITp = 2:  Reading in initial momenta"
        write(6,*)
        write(6,105)"x (au)","y (au)","z (au)"
 105    format(3a12)
        do i=1,natom(im)
          ii=iatom(im)+i
c         read in atomic momenta
          read(5,*)chardum,dum,(pp0(k,ii),k=1,3)
          write(6,106)ii,(pp0(k,ii),k=1,3)
 106      format(i5,4f12.5)
        enddo
        ENDIF
        call MPI_BCAST(ii, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(pp0, 3*ii, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
      elseif (initp(im).eq.-1) then
        IF (my_id.eq.0) THEN
        write(6,*)"INITp = -1:  Initial momenta determined by INITx"
        write(6,*)
        ENDIF
      else
        IF (my_id.eq.0) write(6,*)
     &    "ERROR:  initial conditions flag = ",initp(im),
     &      " is not supported"
        stop
      endif

c INITj
      if (initj(im).eq.0) then
        IF (my_id.eq.0) THEN
        write(6,*)"INITj = 0: Total angular momentum will ",
     &  "not be adjusted"
        write(6,*)
        ENDIF
      elseif (initj(im).eq.1) then
        IF (my_id.eq.0) THEN
        write(6,*)"INITj = 1: Thermal state selected J"
        read(5,*)samptarg(im),letot(im),sampjmin(im),sampjmax(im),
     &   sampjtemp1(im),sampjtemp2(im),sampjbrot1(im),sampjbrot2(im)
        ENDIF
        call MPI_BCAST(samptarg, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
!        call MPI_BCAST(letot, nmol, MPI_LOGICAL,
!     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampjmin, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampjmax, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampjtemp1, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampjtemp2, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampjbrot1, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(sampjbrot2, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
        IF (my_id.eq.0) THEN
        if (letot(im)) 
     &      write(6,*)"Initial total energy is fixed"
        if (.not.letot(im)) 
     &      write(6,*)"Initial Evib+K2 energy is fixed"
        ENDIF
        if (samptarg(im).gt.0.d0) then
          IF (my_id.eq.0) write(6,*)
     &     "Sampled vibrational momenta will be scaled to an ",
     &     "energy of ", samptarg(im)," eV, independent of J"
          samptarg(im)=samptarg(im)/autoev
        elseif (samptarg(im).eq.0.d0) then
          IF (my_id.eq.0) write(6,*)"Sampled momenta will not be scaled"
        elseif (samptarg(im).lt.0.d0) then
          IF (my_id.eq.0) THEN
          write(6,*)"Sampled vibrational momenta will be scaled ",
     &   "according to the rule: Evib(J)=(A+B*J+C*J^2)*D"
          read(5,*)(ejsc(k,im),k=1,4)
          write(6,*)"A=",ejsc(1,im)," eV"
          write(6,*)"B=",ejsc(2,im)," eV"
          write(6,*)"C=",ejsc(3,im)," eV"
          write(6,*)"D=",ejsc(4,im)
          ENDIF
          call MPI_BCAST(ejsc, 4*nmol, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          if (sampjmin(im).lt.0.d0) then
            IF (my_id.eq.0) write(6,*)
     &       "ERROR: Jmin > 0 is required with EJSCALE.f"
            stop
          endif
        endif
        IF (my_id.eq.0) write(6,*)
        if (sampjmin(im).ge.0.d0.and.sampjtemp1(im).lt.0.d0) then
          IF (my_id.eq.0) write(6,*)
     &      "Total angular momentum assigned evenly from ",
     &      sampjmin(im)," to ",sampjmax(im)," hbar"
        elseif (sampjmin(im).ge.0.d0.and.sampjtemp1(im).ge.0.d0) then
          IF (my_id.eq.0) THEN
          write(6,*)"Total angular momentum assigned from a ",
     &     sampjtemp1(im)," to ",sampjtemp2(im),
     &     " K thermal distribution "
          write(6,*)" and from",
     &     sampjmin(im)," to ",sampjmax(im)," hbar with a ",
     &     "rotational constants of ",sampjbrot1(im)," cm-1",
     &     " (x2) and ",sampjbrot2(im)," cm-1 (unique)"
          ENDIF
          sampjbrot1(im)=sampjbrot1(im)/autocmi
          sampjbrot2(im)=sampjbrot2(im)/autocmi
        else
          IF (my_id.eq.0) write(6,*)
     &      "Total angular momentum will not be adjusted"
        endif
        IF (my_id.eq.0) write(6,*)
      else
        IF (my_id.eq.0) write(6,*)
     &     "ERROR:  initial conditions flag INITj= ",initj(im),
     &      " is not supported"
        stop
      endif
c Loop over molecules
      ENDDO

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(nat, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(natom, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(symbol, nat, MPI_CHARACTER,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(mm, nat, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(initx, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(initp, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(initj, nmol, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(ezeroim, nmol, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(natom, nmol, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

c AG orientational information
      IF (my_id.eq.0) THEN
      write(6,*)"AG orientational prescription"
      write(6,*)"-----------------------------"
      ENDIF
      if (nmol.eq.1) then
c        ldofrag=.true.
        if (my_id.eq.0) write(6,*)
     &    "Single AG calculation, no orientational",
     &                   " info required"
        do k=1,3
          comxx(k,1) = 0.d0
          compp(k,1) = 0.d0
        enddo
      else if (nmol.gt.1) then
        IF (my_id.eq.0) THEN
        read(5,*)iorient,ldofrag
        if (ldofrag) then
          write(6,*)"Fragment energies will be computed."
        else
          write(6,*)"Fragment energies will not be computed."
        endif
        write(6,*)
        ENDIF
        if (iorient.eq.0) then
c         read a list of x,y,z and Px,Py,Pz
          do i=1,nmol
            IF (my_id.eq.0) THEN
            read(5,*)(comxx(k,i),k=1,3)
            read(5,*)(compp(k,i),k=1,3)
            write(6,*)"For AG ",i
            write(6,*)"Center of mass coords (x,y,z) (A)  = ",
     &        (comxx(k,i),k=1,3)
            write(6,*)"Center of mass motion (x,y,z) (au) = ",
     &        (compp(k,i),k=1,3)
            ENDIF
            do k=1,3
              comxx(k,i) = comxx(k,i)/autoang
            enddo
            call MPI_BCAST(comxx, 3*nmol, MPI_DOUBLE_PRECISION,
     &                     0, MPI_COMM_WORLD, ierr)      
            call MPI_BCAST(compp, 3*nmol, MPI_DOUBLE_PRECISION,
     &                     0, MPI_COMM_WORLD, ierr)      
          enddo
        elseif (iorient.eq.1) then
          if (nmol.ne.2) then
            IF (my_id.eq.0) write(6,*)
     &       "IORIENT = 1 works with NMOL = 2 only. NMOL = ",
     &       nmol
            stop
          endif
          IF (my_id.eq.0) THEN
          write(6,*)"Quasiclassical bimolecular collisions"
c         quasiclassical state-selected orientation
          read(5,*)rel0qc,tempqc,bminqc,bmaxqc
          write(6,*)"Initial separation (A) = ",rel0qc
          write(6,*)"Impact parameter range (A) = ",bminqc," to ",
     & dabs(bmaxqc)
          if(bmaxqc.gt.0.d0)write(6,*)"Sampled uniformly from 0-bmax"
          if(bmaxqc.lt.0.d0)write(6,*)"Sampled uniformly from 0-bmax^2"
          ENDIF
          call MPI_BCAST(rel0qc, 1, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(tempqc, 1, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(bminqc, 1, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          call MPI_BCAST(bmaxqc, 1, MPI_DOUBLE_PRECISION,
     &                   0, MPI_COMM_WORLD, ierr)      
          if (tempqc.gt.0.d0) then
            IF (my_id.eq.0) write(6,*)
     &      "Relative temperature (K) = ",tempqc
          else if (tempqc.eq.0.d0) then
            IF (my_id.eq.0) write(6,*)
     &      "Relative temperature set to sampled ",
     &      "J temperature for AG 1"
          else if (tempqc.lt.0.d0) then
            IF (my_id.eq.0) write(6,*)
     &        "Constant initial relative energy = ",
     &                 dabs(tempqc)," eV"
            tempqc=tempqc/autoev
          endif
          if (dabs(bmaxqc).gt.rel0qc) then
             IF (my_id.eq.0) write(6,*)
     &          "Maximum impact parameter (",bmaxqc," A) cannot",
     &          " be greater than initial separation (",rel0qc," A)"
             stop
          endif
          bminqc=bminqc/autoang
          bmaxqc=bmaxqc/autoang
          rel0qc=rel0qc/autoang
        else
          IF (my_id.eq.0) write(6,*)
     &      "IORIENT = ",iorient," is not supported"
          stop
        endif
      else
        IF (my_id.eq.0) write(6,*)"Cant have AG < 0"
        stop
      endif
      IF (my_id.eq.0) write(6,*)

c read termination condition
      IF (my_id.eq.0) THEN
      write(6,*)"Termination condition"
      write(6,*)"---------------------"
      read(5,*)termflag,t_nstep
c      read(5,*)termflag,t_nstep,lwell
      ENDIF
      call MPI_BCAST(termflag, 1, MPI_INTEGER,
     &               0, MPI_COMM_WORLD, ierr)      
      call MPI_BCAST(t_nstep, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)      
      lwell=0
      if (termflag.eq.0) then
c       fixed number of steps per trajectory
        IF (my_id.eq.0) write(6,*)
     &   "TERMFLAG = 0:  Fixed number of steps per trajectory"
      elseif (termflag.eq.1) then
c       fixed simulation time
        IF (my_id.eq.0) THEN
        read(5,*)t_stime
        write(6,*)"TERMFLAG = 1:  Trajectories will run for ",
     &  t_stime," fs"
        t_stime=t_stime/autofs
        ENDIF
        call MPI_BCAST(t_stime, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
      elseif (termflag.eq.2) then
c       converge gradient
        IF (my_id.eq.0) THEN
        read(5,*)t_gradmag
        write(6,*)"TERMFALG = 2:  Trajectories will terminate ",
     &  "when the magnitude of the gradient"
        write(6,*)"               is less than ",
     &  t_gradmag," eV/A"
        t_gradmag=t_gradmag*autoang/autoev
        ENDIF
        call MPI_BCAST(t_gradmag, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
      elseif (termflag.eq.3) then
c       dissociation
        IF (my_id.eq.0) THEN 
        write(6,*)"TERMFLAG = 3:  Association/dissociation",
     &   " termination conditions"
        write(6,*)"Outcome label      Bond        Distance (A)"
        read(5,*)t_noutcome
        do i=1,t_noutcome
          read(5,*)t_symb(i,1),t_symb(i,2),t_r(i)
          if (t_r(i).ge.0.d0) then
          write(6,200)i,t_symb(i,1),t_symb(i,2),t_r(i),"   dissociation"
          else
          write(6,200)i,t_symb(i,1),t_symb(i,2),-t_r(i),"   association"
          endif
 200      format(i5,15x,a2,"-",a2,f15.5,a)
          t_r(i)=t_r(i)/autoang
        enddo
        write(6,*)
        ENDIF
        call MPI_BCAST(t_noutcome, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)      
        call MPI_BCAST(t_r, t_noutcome, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)      
      elseif (termflag.eq.4) then
c     monitor dmag
        IF (my_id.eq.0) write(6,*)
     &    "Terminating trajectory at the first minimum in |d|."
      else
        IF (my_id.eq.0) write(6,*)
     &    "Dont know how to handle TERMFLAG = ",termflag
        stop
      endif

      IF (my_id.eq.0) THEN
      if (t_nstep.gt.0) then
        write(6,*)"The maximum number of steps per trajectory is ",
     &  t_nstep
      else
        write(6,*)"The maximum number of steps per trajectory is ",
     &  "unlimited"
      endif
      write(6,*)
      ENDIF

      if (lwell.gt.0) write(6,*)"Selecting geoms in PE well ",lwell

c output flags
      IF (my_id.eq.0) THEN
      write(6,*)"Output options"
      write(6,*)"--------------"
      read(5,*)ii,(itmp(k),k=1,ii)
      if (ii.eq.0) then 
      write(6,*)"Writing to all output files"
      do i=1,99
        lwrite(i) = .true.
      enddo
      else
      write(6,*)"Writing to units ",(itmp(k),k=1,ii)
      do i=1,99
        lwrite(i) = .false.
      enddo
      do i=1,ii
        lwrite(itmp(i)) = .true.
        if (itmp(i).eq.77) 
     & write(6,*)"Unit 77 chosen: Writing special named output files"
      enddo
      endif
      write(6,*)
      ENDIF

c print headers
      IF (my_id.eq.0) THEN
      if (lwrite(80)) then
        write(80,*)"[Molden Format]"
        write(80,*)"[GEOMETRIES] XYZ"
      endif
      if (lwrite(81)) then
        write(81,*)"title"
        write(81,*)"momenta XYZ au "
      endif
      nrecl=(3+5*nat)*10
      if (lwrite(82)) open(82,form='unformatted',
     &   status='unknown',access='direct',RECL=nrecl)
      if (lwrite(83)) open(83,form='unformatted',
     &   status='unknown',access='direct',RECL=nrecl)
      
      do im=1,nmol
      if (lems(im).or.lwell.gt.0) then
        if (lwell.gt.0) then
          nrecl=(3+5*nat)*10
        else
          nrecl=(3+5*natom(im))*10
        endif
        if (lwrite(86)) then
          write(86,*)"[Molden Format]"
          write(86,*)"[GEOMETRIES] XYZ"
        endif
        if (lwrite(87)) then
          write(87,*)"title"
          write(87,*)"momenta XYZ au "
        endif
        if (lwrite(88)) open(88,form='unformatted',
     &   status='unknown',access='direct',RECL=nrecl)
        if (lwrite(89)) open(89,form='unformatted',
     &     status='unknown',access='direct',RECL=nrecl)
      endif
      enddo
      ENDIF

c check for bad input combinations
      if (potflag.eq.1.and.nat.ne.3) then
        IF (my_id.eq.0)write(6,*)"POTFLAG = 1 must be used with 3 atoms"
        stop 
      endif

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      return

      end
