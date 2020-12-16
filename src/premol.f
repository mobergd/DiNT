      subroutine premol(im)

c Precompute information for each atom group.

      implicit none
      include 'param.f'
      include 'c_sys.f'

c MPI
      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)
 
      double precision xxm(3,mnat),ppm(3,mnat),mmm(mnat)
      double precision nmvecm(3,mnat,3*mnat),freqm(3*mnat),
     & rturnm(3*mnat),ewellm,nmqnm(3*mnat)
      integer im,i,j,ii,k
      character*2 symb(mnat)

ccccc MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

      print *,"Got to premol.f on proc ",my_id

c MPI
      IF (my_id.eq.0) THEN
        write(6,*)"Precomputing info for AG ",im
        write(6,*)"------------------------------"
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

      if (initx(im).eq.2.or.initx(im).eq.5) then
c       SET UP FOR NORMAL MODE METHODS
c       set to initial structure
        do i=1,natom(im)
          ii=i+iatom(im)
          mmm(i) = mm(ii)
          do j=1,3
            xxm(j,i) = xx0(j,ii)
          enddo
        enddo

        do i=1,natom(im)*3
          nmqnm(i) = nmqn(i,im)
        enddo
        call normod(xxm,mmm,natom(im),nmvecm,freqm,ewellm,
     &   repflag,nsurf0,rturnm,nmqnm,lreadhess,ldofrag,
     &   symb,nmtype(im))
        ewell(im) = ewellm
        do i=1,3*natom(im)
          freq(i,im)=freqm(i)
          rturn(i,im)=rturnm(i)
          do j=1,3
            do k=1,natom(im)
              nmvec(j,k,i,im)=nmvecm(j,k,i)
            enddo
          enddo
        enddo
        if (sampwell(im).eq.2) then
c       NORMAL MODE METHOD for second PE well
c       set to initial structure
          do i=1,natom(im)
            ii=i+iatom(im)
            mmm(i) = mm(ii)
            do j=1,3
              xxm(j,i) = xx02(j,ii)
            enddo
          enddo

          do i=1,natom(im)*3
            nmqnm(i) = nmqn(i,im)
          enddo
          call normod(xxm,mmm,natom(im),nmvecm,freqm,ewellm,
     &     repflag,nsurf0,rturnm,nmqnm,lreadhess,ldofrag,
     &     symb,nmtype(im))
          ewell2(im) = ewellm
          do i=1,3*natom(im)
            freq2(i,im)=freqm(i)
            rturn2(i,im)=rturnm(i)
            do j=1,3
              do k=1,natom(im)
                nmvec2(j,k,i,im)=nmvecm(j,k,i)
              enddo
            enddo
          enddo
        endif
      else if (initx(im).eq.3) then
c ATOM-DIATOM METHOD
        IF (my_id.eq.0) THEN
        write(6,*)"INITx = 3:  Atom-diatom initial conditions"
        ENDIF
        call preatomdiatom(escatad,jjad,vvad,arrad,
     &    rrad,mmm,rinad,routad,tauad,bmaxad,
     &    bminad,pprelad,ecolad,nsurf0,nsurft,easym,rasym)
      else if (initx(im).eq.4) then
c DIATOM METHOD
        IF (my_id.eq.0) THEN
        write(6,*)"INITx = 4:  Diatom initial conditions"
        ENDIF
        call prediatom(vvdi(im),jjdi(im),rmindi(im),mmm,nsurf0,
     &     symb,taudi(im),rindi(im),routdi(im))
      else
        IF (my_id.eq.0) THEN
        write(6,*)"Nothing to do..."
        write(6,*)
        ENDIF
      endif

      return

      end
