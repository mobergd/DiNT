      subroutine initmol(im)

c Initialize AG #IM.  Use information computed and stored from PREMOL
c to set up each trajectory.

      implicit none
      include 'param.f'
      include 'c_traj.f'
      include 'c_sys.f'
 
c MPI
      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      double precision xtot,ytot,ztot,mtot,timpnow,kinow,etmp
!      double precision xxm1(mnat),xxm2(mnat),xxm3(mnat),
!     &  ppm1(mnat),ppm2(mnat),ppm3(mnat),
!     &  xx1(mnat),xx2(mnat),xx3(mnat),pp1(mnat),pp2(mnat),pp3(mnat)
      double precision xxm(3,mnat),ppm(3,mnat),mmm(mnat),
     &  phop(mnsurf),dmag,ppj(3,mnat),ek2,ejrot
      integer im,i,j,k,ii
      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & pema(mnsurf),gpema(3,mnat,mnsurf),dvec(3,mnat,mnsurf,mnsurf)
      character*2 symb(mnat)

c     local
      integer nqn,whwell,ntry
      double precision molpe,etarget,scalef,sampjj,sampkk,eig(3),etest,
     & dum3(3)
c      logical letot

c HARD-CODED --- changed to input flag letot(im)
c      letot = .true. ! fixed total energy
c      letot = .false. ! fixed vib+K2 energy

ccccc MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

      IF (my_id.eq.0) THEN
      write(6,*)"Generating initial conditions for AG ",im
      write(6,*)"-------------------------------------------"
      ENDIF

c     store coords for this molecule in their own array for easier manipulation
      do i=1,natom(im)
        ii=i+iatom(im)
        mmm(i) = mm(ii)
        symb(i) = symbol(ii)
        do j=1,3
          xxm(j,i) = 0.d0
          ppm(j,i) = 0.d0
        enddo
      enddo
      molpe=0.d0
      ntry=0

c GENERATE INITIAL COORDINATES
 111  continue
      if (initx(im).eq.0) then
c       they were read in readin and stored in xx0
        IF (my_id.eq.0)
     &    write(6,*)"INITx = 0:  Initial coordinates are ",
     &      "specified by the user"
        do i=1,natom(im)
        do j=1,3
          ii=i+iatom(im)
          xxm(j,i) = xx0(j,ii)
        enddo
        enddo
      elseif (initx(im).eq.1) then
        IF (my_id.eq.0)
     &    write(6,*)"INITx = 1:  Initial coordinates are ",
     &    "selected randomly"
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call ranclu(xxm,natom(im))
      elseif (initx(im).eq.2) then
        IF (my_id.eq.0)
     &    write(6,*)"INITx = ",initx(im),":  Normal mode analysis"
c       assign minE structure
        do i=1,natom(im)
        do j=1,3
          ii=i+iatom(im)
          xxm(j,i) = xx0(j,ii)
          ppm(j,i) = 0.d0
        enddo
        enddo
c       populate normal modes according to input quantum numbers
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call popnorm(xxm,ppm,mmm,natom(im),im)
      elseif(initx(im).eq.5.and..not.lems(im)) then
c       pick normal mode quanta randomly & assign initial structure
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call vibwells(xxm,ppm,scale0im(im),im)
c       populate normal modes according to the quantum numbers from vibwells
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call popnorm(xxm,ppm,mmm,natom(im),im)
      elseif (initx(im).eq.5.and.lems(im)) then 
        if (itraj.gt.1) then
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call hessems(xxm,ppm,molpe,mmm,symb,nbrea(im),im) 
        elseif (itraj.eq.1) then 
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call hessems(xxm,ppm,molpe,mmm,symb,ninc(im),im)
        endif
      elseif (initx(im).eq.3) then
        IF (my_id.eq.0)
     &    write(6,*)"INITx = 3:  Atom-diatom initial conditions "
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call atomdiatom(arrad,rrad,ecolad,bminad,
     &   bmaxad,jjad,rinad,routad,tauad,
     &   pprelad,mmm,xxm,ppm)
      elseif (initx(im).eq.4) then
        IF (my_id.eq.0)
     &    write(6,*)"INITx = 4:  Diatom initial conditions "
        call diatom(jjdi(im),taudi(im),rindi(im),routdi(im),mmm,xxm,
     &       ppm)
      elseif (initx(im).eq.6.and..not.lems(im)) then
        IF (my_id.eq.0)
     &    write(6,*)"INITx = 6:  Read from a separate file"
        if (sampwell(im).eq.2) then
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call ranwell(whwell,relwell1(im))
          if (whwell.eq.1) then
            call ransamp(lbinsamp(im),lems(im),samptot(im),
     &       sampfilexx(im),sampfilepp(im),natom(im),xxm,ppm,molpe)
          else if (whwell.eq.2) then
            call ransamp(lbinsamp2(im),lems(im),samptot2(im),
     &       sampfilexx2(im),sampfilepp2(im),natom(im),xxm,ppm,molpe)
          endif
        else 
          call ransamp(lbinsamp(im),lems(im),samptot(im),
     &     sampfilexx(im),sampfilepp(im),natom(im),xxm,ppm,molpe)
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call MPI_BCAST(molpe, 1, MPI_DOUBLE_PRECISION,
     &                 0, MPI_COMM_WORLD, ierr)
          call MPI_BCAST(ppm, 3*natom(im), MPI_DOUBLE_PRECISION,
     &             0, MPI_COMM_WORLD, ierr)
        endif
      elseif (initx(im).eq.6.and.lems(im)) then
        if (itraj.gt.1) then
          call ems(xxm,ppm,molpe,mmm,symb,nbrea(im),im)
        elseif (itraj.eq.1) then          
          call ems(xxm,ppm,molpe,mmm,symb,ninc(im),im)
        endif
      else
        IF (my_id.eq.0)
     &    write(6,*)"ERROR:  INITx = ",initx(im)," is not an option"
        stop
      endif

c put center of mass at the origin
      IF (my_id.eq.0) THEN
      write(6,100)"Placing CoM at origin"
      write(6,*)
      ENDIF
 100  format(1x,a,3f12.5,1x,a)
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,natom(im)
        xtot = xtot + mmm(i)*xxm(1,i)
        ytot = ytot + mmm(i)*xxm(2,i)
        ztot = ztot + mmm(i)*xxm(3,i)
        mtot = mtot + mmm(i)
      enddo
      do i=1,natom(im)
        xxm(1,i) = xxm(1,i) - xtot/mtot
        xxm(2,i) = xxm(2,i) - ytot/mtot
        xxm(3,i) = xxm(3,i) - ztot/mtot
      enddo

c spin randomly if necessary
c HACK AJ
      if (iorient.eq.1) call spin(xxm,ppm,natom(im))
      if (iorient.eq.2) call spin(xxm,ppm,natom(im))
c HACK AJ

c print coordinates
c      if (.false.) then ! Minimal printing
c MPI
      IF (my_id.eq.0) THEN
      write(6,*)"Initial coordinates"
      write(6,105)"    # symb","x (A)","y (A)","z (A)"
      ENDIF
 105  format(a10,3a12)
      do i=1,natom(im)
        ii=iatom(im)+i
        IF (my_id.eq.0) 
     &   write(6,101)ii,symbol(ii),xxm(1,i)*autoang,
     &            xxm(2,i)*autoang,xxm(3,i)*autoang
 101    format(i5,a5,3f12.5)
      enddo
      IF (my_id.eq.0) write(6,*)
c      endif

c if we haven't calcualted the PE yet 
c      if (molpe.eq.0.d0) then
      call MPI_BCAST(ldofrag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      if (ldofrag) then
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call getpem(xxm,natom(im),pema,pemd,gpema,gpemd,dvec,symb)
        if (repflag.eq.1) then
          molpe = pemd(nsurf,nsurf)
        else
          molpe = pema(nsurf)
        endif
c       correct for this atom groups zero of energy
        molpe = molpe + ezero - ezeroim(im)
      else
        IF (my_id.eq.0)
     &    write(6,*)"Initial potential energy evalulation skipped"
        molpe = 0.d0
      endif
c      endif

c GENERATE INITIAL MOMENTA
      if (initp(im).eq.0) then
c     random thermal sample
        IF (my_id.eq.0)
     &   write(6,*)"INITp = 0:  Generating initial velocities from a",
     &   " random thermal distribution"
        call rantherm(ppm,mmm,natom(im),temp0im(im))
        IF (my_id.eq.0) 
     &    write(6,106)"Desired temp    = ",temp0im(im)," K"
      elseif (initp(im).eq.1) then
        IF (my_id.eq.0)
     &    write(6,*)"INITp = 1:  Setting initial momenta to zero"
        do i=1,natom(im)
          do j=1,3
            ppm(j,i)=0.d0
          enddo
        enddo
      elseif (initp(im).eq.2) then
        IF (my_id.eq.0)write(6,*)
     &   "INITp = 2:  Initial momenta are specified by the user"
        do i=1,natom(im)
          do j=1,3
            ppm(j,i)=pp0(j,i)
          enddo
        enddo
      elseif (initp(im).eq.-1) then
        IF (my_id.eq.0) THEN
        write(6,*)"INITp = -1:  Inititial momenta selected along with",
     & " initial coordinates"
        write(6,*)
        ENDIF
      else
        IF (my_id.eq.0)
     &    write(6,*)"ERROR:  INITp = ",initp(im)," is not an option"
        stop
 106  format(1x,a,f12.5,1x,a)
      IF (my_id.eq.0) write(6,*)
      endif

c BCAST pp matrix
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ppm, 3*natom(im), MPI_DOUBLE_PRECISION,
     &             0, MPI_COMM_WORLD, ierr)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call gettemp(ppm,mmm,natom(im),temp,ke)
      IF (my_id.eq.0) THEN
      write(6,106)"Calculated temp = ",temp," K"
      write(6,106)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF


c MANIPULATE COORDINATES AND MOMENTA
c     remove overall momenta
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,natom(im)
        xtot=xtot+ppm(1,i)
        ytot=ytot+ppm(2,i)
        ztot=ztot+ppm(3,i)
        mtot = mtot + mmm(i)
      enddo
      do i=1,natom(im)
        ppm(1,i)=ppm(1,i)-xtot*mmm(i)/mtot
        ppm(2,i)=ppm(2,i)-ytot*mmm(i)/mtot
        ppm(3,i)=ppm(3,i)-ztot*mmm(i)/mtot
      enddo
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call gettemp(ppm,mmm,natom(im),temp,ke)
      IF (my_id.eq.0) THEN
      write(6,*)"Removed CoM motion"
      write(6,106)"Calculated temp = ",temp," K"
      write(6,106)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF

      if (initj(im).eq.1) then

C       SCALE MOMENTUM
        if (sampjmin(im).ge.0.d0) then
          IF (my_id.eq.0) print *,"Removing angular momentum"
          call MPI_BARRIER(MPI_COMM_WORLD, ierr)
          call noang(xxm,ppm,mmm,natom(im))
          call gettemp(ppm,mmm,natom(im),temp,ke)
          IF (my_id.eq.0) THEN
          write(6,106)"Calculated temp = ",temp," K"
          write(6,106)"Total KE        = ",ke*autoev," eV"
          write(6,*)
          ENDIF
        endif
 
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call gettemp(ppm,mmm,natom(im),temp,ke)
c         print *,"Unscaled total energy = ",(molpe+ke)*autoev," eV"
c         print *

        if (samptarg(im).gt.0.d0) then
          etarget=samptarg(im)  ! fixed-energy scaling
c          print *,"Target total energy = ",etarget*autoev," eV"
c          print *
        elseif (samptarg(im).le.0.d0) then  ! random J
          IF (my_id.eq.0) write(6,*)"Choosing J"
          if (sampjtemp(im).ge.0.d0) then
            call ranno(sampjtemp(im),sampjtemp1(im),sampjtemp2(im))
            IF (my_id.eq.0) write(6,*)
     &        "Selecting J from distribution at ",sampjtemp(im),"K"
          endif
          call ranj(sampjmin(im),
     &      sampjmax(im),sampjtemp(im),sampjj,sampkk,
     &      sampjbrot1(im),sampjbrot2(im))
          if (samptarg(im).lt.0.d0) then
            call ejscale(im,sampjj,ejsc,etarget)
          else
            etarget=0.d0
          endif
c          print *,"Target total energy = ",etarget*autoev," eV"
c          print *
        endif

        ek2=0.d0
        do i=1,natom(im)
          do j=1,3
            ppj(j,i)=0.d0
          enddo
        enddo
        if (sampjmin(im).ge.0.d0) then
c          write(6,*)"Computing rotational energy and momenta"
          call addrot(ppj,xxm,mmm,natom(im),sampjj,sampkk,
     &         sampjbrot1(im),sampjbrot2(im),ek2,ejrot,dum3)
        endif

        IF (my_id.eq.0) THEN
        write(6,*)"Initial kinetic energy          = ",ke*autoev
        write(6,*)"Initial potential energy        = ",molpe*autoev
        write(6,*)"Initial K2 energy               = ",ek2*autoev
        write(6,*)"Initial rot energy              = ",ejrot*autoev
        write(6,*)"Target energy                   = ",etarget*autoev
        ENDIF
        if (letot(im)) etest=ejrot
        if (.not.letot(im)) etest=ek2
        if ((molpe+etest).gt.etarget.and.samptarg(im).ne.0.d0) then
c          if (.false.) then
          IF (my_id.eq.0) THEN
          print *,"Potential energy + Trot(",(molpe+etest)*autoev,
     &       " eV) > Target total energy (",etarget*autoev,") eV"
          print *,"Picking another geometry"
          ENDIF
          molpe=0.d0
c          stop ! TEMP AJ
          ntry=ntry+1
          if (ntry.eq.10) stop
          go to 111
        endif

        if (samptarg(im).ne.0.d0) then
          scalef=dsqrt((etarget-molpe-etest)/ke)
c          scalef=dsqrt((etarget-etest)/ke)
        else
          scalef=1.d0
        endif
        IF (my_id.eq.0) THEN
        print *,"Scaling factor",scalef
        print *
        ENDIF
        do i=1,natom(im)
          do j=1,3
            ppm(j,i)=ppm(j,i)*scalef
          enddo
        enddo

        call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        call gettemp(ppm,mmm,natom(im),temp,ke)
        IF (my_id.eq.0) THEN
        write(6,*)"After scaling P"
        write(6,106)"Calculated temp = ",temp," K"
        write(6,106)"Total KE        = ",ke*autoev," eV"
        write(6,106)"Total KE+V      = ",(molpe+ke)*autoev," eV"
        write(6,106)"Total KE+V+EK2  = ",(molpe+ke+ek2)*autoev," eV"
        write(6,*)
        write(6,*)"Adding rotation"
        ENDIF
        do i=1,natom(im)
          do j=1,3
            ppm(j,i)=ppm(j,i)+ppj(j,i)
          enddo
        enddo
      endif

c     calculate angular motion
      if (natom(im).gt.2.and.initx(im).ne.3) then
c     don't do for diatoms or for special atom-diatom initial conditions
        call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        call ange(xxm,ppm,mmm,natom(im),eig,bigj,bigjtot,erot,erottot)
        IF (my_id.eq.0) THEN
        write(6,106)"Total angular momentum      ",bigjtot," au"
        write(6,107)"Angular momentum components ",bigj," au"
        write(6,106)"Total rotational energy     ",erottot*autoev," eV"
        write(6,107)"Rotational energy components",erot(1)*autoev,
     &   erot(2)*autoev,erot(3)*autoev," eV"
        write(6,*)
        ENDIF
 107  format(1x,a,3f12.5,1x,a)

      if (initj(im).ne.1) then
c Comment \/ this to skip removal of angular momentum
c      call noang(xxm,ppm,mmm,natom(im))
c      write(6,*)"Removed overall angular motion"
c      call ange(xxm,ppm,mmm,natom(im),eig,bigj,bigjtot,erot,erottot)
c      call gettemp(ppm,mmm,natom(im),temp,ke)
c      write(6,106)"Total angular momentum      ",bigjtot," au"
c      write(6,107)"Angular momentum components ",bigj," au"
c      write(6,106)"Total rotational energy     ",erottot*autoev," eV"
c      write(6,107)"Rotational energy components",erot(1)*autoev,
c     &   erot(2)*autoev,erot(3)*autoev," eV"
c Comment ^ this to skip removal of angular momentum
        call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        call gettemp(ppm,mmm,natom(im),temp,ke)
        IF (my_id.eq.0) THEN
        write(6,*)
        write(6,106)"Calculated temp = ",temp," K"
        write(6,106)"Total KE        = ",ke*autoev," eV"
        write(6,*)
        ENDIF
      endif

      if (initp(im).eq.0.and.escale0im(im).gt.0.d0) then
c       rescale momenta
        do i=1,natom(im)
          do j=1,3
            ppm(j,i)=ppm(j,i)*dsqrt(escale0im(im)/ke)
          enddo
        enddo
        call MPI_BARRIER(MPI_COMM_WORLD, ierr) 
        call gettemp(ppm,mmm,natom(im),temp,ke)
        IF (my_id.eq.0) THEN
        write(6,*)"Rescaling KE to ",escale0im(im)*autoev," eV"
        write(6,106)"Calculated temp = ",temp," K"
        write(6,106)"Total KE        = ",ke*autoev," eV"
        write(6,*)
        ENDIF
      endif

      if (.false.) then ! Minimal printing
        IF (my_id.eq.0) THEN
        write(6,*)"Initial momenta for this AG"
        write(6,105)"    # symb","x (au)","y (au)","z (au)"
        do i=1,natom(im)
          ii=i+iatom(im)
          write(6,102)ii,symbol(ii),ppm(1,i),ppm(2,i),ppm(3,i)
        enddo
        write(6,*)
        ENDIF
      endif

 102  format(i6,a5,3e12.5)

      else

      bigjtot = 0.d0
      erottot = 0.d0

      endif

c print coordinates
      if (.false.) then ! Minimal printing
        IF (my_id.eq.0) THEN
        write(6,*)"Initial coordinates"
        write(6,105)"    # symb","x (A)","y (A)","z (A)"
        do i=1,natom(im)
          ii=iatom(im)+i
          write(6,101)ii,symbol(ii),xxm(1,i)*autoang,
     &             xxm(2,i)*autoang,xxm(3,i)*autoang
        enddo
        write(6,*)
        ENDIF
      endif

c WRITE INITIAL ENERGY AND ANGULAR MOMENTUM
      IF (my_id.eq.0) write(6,*)"Info for this AG only"
      call MPI_BCAST(ldofrag, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      if (ldofrag) then
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call getpem(xxm,natom(im),pema,pemd,gpema,gpemd,dvec,symb)
        if (repflag.eq.1) then
          pe = pemd(nsurf,nsurf)
        else
          pe = pema(nsurf)
        endif
c       correct for this atom groups zero of energy
        pe = pe + ezero - ezeroim(im)
      else
        IF (my_id.eq.0) write(6,*)
     &    "Initial potential energy evalulation skipped"
        pe = 0.d0
      endif
      peagi(im)=pe
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call gettemp(ppm,mmm,natom(im),temp,ke)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      IF (my_id.eq.0) THEN
      call ange(xxm,ppm,mmm,natom(im),eig,bigj,bigjtot,erot,erottot)
      ENDIF
      call MPI_BCAST(bigj, 3, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(bigjtot, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(erot, 3, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(erottot, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      te = pe+ke
      IF (my_id.eq.0) THEN   
      write(6,106)"Initial potential energy   = ",
     &   pe*autoev," eV"
      write(6,106)"Initial kinetic energy     = ",
     &   ke*autoev," eV"
      write(6,106)"Initial rotational energy  = ",
     &   erottot*autoev," eV"
      write(6,106)"Initial vibrational energy = ",
     &   (ke-erottot)*autoev," eV"
      write(6,106)"Initial total energy       = ",
     &   te*autoev," eV"
      write(6,106)"Initial temp               = ",
     &   temp," K"
      write(6,106)"Initial angular momentum   = ",
     &   bigjtot," au"
      write(6,*)
      ENDIF

c save AG info
      erottotagi(im)=erottot
      evibagi(im)=ke-erottot
      eintagi(im)=ke
      bigjtotagi(im)=bigjtot
      do i=1,3
        bigjagi(i,im)=bigj(i)
        erotagi(i,im)=erot(i)
      enddo

c write manipulated data back to xx and pp arrays
      IF (my_id.eq.0) THEN
      do i=1,natom(im)
        ii=i+iatom(im)
        do j=1,3
          xx(j,ii) = xxm(j,i)
          pp(j,ii) = ppm(j,i)
        enddo
      enddo
      ENDIF

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(xx, 3*natom(im), MPI_DOUBLE_PRECISION,
     &             0, MPI_COMM_WORLD, ierr)
      call MPI_BCAST(pp, 3*natom(im), MPI_DOUBLE_PRECISION,
     &             0, MPI_COMM_WORLD, ierr)

c temp
      rtrans=0.d0

      call MPI_BARRIER(MPI_COMM_WORLD, ierr) 

      return

      end
