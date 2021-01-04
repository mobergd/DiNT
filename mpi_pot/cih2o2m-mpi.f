      subroutine pot(symb,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

c NEW "AUTOFIT" PESs

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension x2(maxatom),y2(maxatom),z2(maxatom)
      dimension x3(maxatom),y3(maxatom),z3(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      dimension dvdx2(maxatom),dvdy2(maxatom),dvdz2(maxatom)
      dimension dvdx3(maxatom),dvdy3(maxatom),dvdz3(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      character*2 symb(maxatom),symb2(maxatom),symb3(maxatom)
      integer at(maxatom),at2(maxatom),at3(maxatom)
      integer itinkxyz
      logical ltinkm

      dimension tmpprint(50)
      common/tmp/tmpprint

ccccc MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

      v=0.d0
      v2=0.d0
      v3=0.d0
      do i=1,natom
        symb2(i)="xx"
        symb3(i)="xx"
        x2(i)=0.d0
        y2(i)=0.d0
        z2(i)=0.d0
        x3(i)=0.d0
        y3(i)=0.d0
        z3(i)=0.d0
        dvdz2(i)=0.d0
        dvdy2(i)=0.d0
        dvdx2(i)=0.d0
        dvdz3(i)=0.d0
        dvdy3(i)=0.d0
        dvdx3(i)=0.d0
        dvdz(i)=0.d0
        dvdy(i)=0.d0
        dvdx(i)=0.d0
      enddo

      do i=1,natom
        at(i)=0
        if ((symb(i).eq."H").or.
     &    (symb(i).eq."h"))  at(i)=1       ! hydrogen
        if ((symb(i).eq."C").or.
     &    (symb(i).eq."c"))  at(i)=2       ! carbon
        if ((symb(i).eq."O").or.
     &    (symb(i).eq."o"))  at(i)=3       ! o
c        if ((symb(i).eq."O1").or.
c     &    (symb(i).eq."o1"))  at(i)=4       ! h2o2
c        if ((symb(i).eq."H1").or.
c     &    (symb(i).eq."H1"))  at(i)=5       ! h2o2
        if ((symb(i).eq."He").or.
     &    (symb(i).eq."he").or.
     &    (symb(i).eq."HE")) at(i)=21      ! helium
        if ((symb(i).eq."Ar").or.
     &    (symb(i).eq."ar").or.
     &    (symb(i).eq."AR")) at(i)=23      ! argon
        if ((symb(i).eq."O1").or.
     &    (symb(i).eq."o1")) at(i)=11      ! 
        if ((symb(i).eq."H1").or.
     &    (symb(i).eq."h1")) at(i)=12      ! 
        if (at(i).eq.0) then ! atom not found
           IF (my_id.eq.0)
     &     write(6,*)"Atom # ",i," (",symb(i),") not found"
           stop
        endif
      enddo

      natom2=0
      natom3=0
      natom4=0
      ii=0
      ltinkm=.false.
      do i=1,natom
        if (at(i).le.10) then ! ch2oo
          ii=ii+1
          natom2=natom2+1
          x2(ii)=x(i)
          y2(ii)=y(i)
          z2(ii)=z(i)
          symb2(ii)=symb(i)
          at2(ii)=at(i)
        elseif (at(i).le.20) then ! h2o2
          ii=ii+1
          natom3=natom3+1
          x2(ii)=x(i)
          y2(ii)=y(i)
          z2(ii)=z(i)
          symb2(ii)=symb(i)
          at2(ii)=at(i)
        else
          natom4=natom4+1 ! m
        endif
      enddo

      if (natom3.ne.0.and.natom2.ne.0) 
     &  call lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,ii,maxatom) !  Interaction PES  ch2oo+h2o2
!      print *,"Step 1 v = ",v,"on proc",my_id
!      print *,"Step 1 v2 = ",v2,"on proc",my_id
!      print *,"Step 1 v3 = ",v3,"on proc",my_id
      if (natom3.ne.0.and.natom2.ne.0.and.natom4.ne.0) 
     &  call rgexp(at,x,y,z,v3,dvdx3,dvdy3,dvdz3,natom,maxatom) !  Interaction PES  ch2oo+h2o2 + m
!      print *,"Step 2 v = ",v,"on proc",my_id
!      print *,"Step 2 v2 = ",v2,"on proc",my_id
!      print *,"Step 2 v3 = ",v3,"on proc",my_id

      itinkxyz=0
      if (natom3.eq.0) itinkxyz=1   ! only targets
      if (natom2.eq.0) itinkxyz=2   ! only baths
      if (ii.ne.0)
     &    call tinkerpot(symb2,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,
     &          ii,maxatom,itinkxyz)
c      if (.not.ltinkm.and.natom3.ne.0)
c     &   call bath(at3,x3,y3,z3,v3,dvdx3,dvdy3,dvdz3,natom3,maxatom) !  1D diatomic baths

!      print *,"Step 3 v = ",v,"on proc",my_id
!      print *,"Step 3 v2 = ",v2,"on proc",my_id
!      print *,"Step 3 v3 = ",v3,"on proc",my_id

      tmpprint(1)=v3   ! interaction
      tmpprint(3)=v2   ! target internal
      tmpprint(4)=0.   ! bath internal
      tmpprint(11)=v    ! interaction

      v=v+v2+v3

      do i=1,natom
        dvdx(i)=dvdx(i)+dvdx2(i)+dvdx3(i)
        dvdy(i)=dvdy(i)+dvdy2(i)+dvdy3(i)
        dvdz(i)=dvdz(i)+dvdz2(i)+dvdz3(i)
      enddo

c hack int PES
      if (natom3.ne.0.and.natom2.ne.0) then
        do i=1,ii
          x2(i)=x(i)
          y2(i)=y(i)
          z2(i)=z(i)
          if (i.ge.6) x2(i)=x2(i)+20.
        enddo
        call lsfit(at,x2,y2,z2,v2,dvdx2,dvdy2,dvdz2,ii,maxatom) !  Interaction PES
        tmpprint(11)=tmpprint(11)-v2    ! interaction
      endif

      return

      end

! ONE DIMENSIONAL DIATOMIC BATHS
      subroutine bath(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)
      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      integer at(maxatom)
      parameter(autocmi=219474.63067d0)
      parameter(autoang=0.529177249d0)
      parameter(autoev=27.2113961d0)

      v=0.
      do i=1,natom
        dvdx(i)=0.d0
        dvdy(i)=0.d0
        dvdz(i)=0.d0
      enddo

      if (natom.eq.1) return

      if (natom.gt.2) then
        print *,"Can't handle more than 2 bath atoms"
        stop
      endif

      if (natom.eq.2) then
        dx=x(1)-x(2)
        dy=y(1)-y(2)
        dz=z(1)-z(2)
        rr=dsqrt(dx*dx+dy*dy+dz*dz)

        if (at(1).eq.25.and.at(2).eq.25) then
! H2 bath
!       From Hack's fit (eq 8 in Hack, Truhlar, JCP 110, 4315 (1999))
!                        to Kolos and Wolniewicz JCP 43, 2429 (1965)
!       Rmin = 1.40121 au, Vmin = -4.74772265 eV relative to H+H
        c1=139.7160d0        ! eV
        c2=-123.8978d0       ! eV / bohr
        c3=3.4031d0          ! 1 / bohr
        c4=-6.8725d0         ! eV / bohr**2
        c5=-23.0440d0        ! eV / bohr**3
        c6=2.032d0           ! 1 / bohr

        v=(c1+c2*rr)*dexp(-c3*rr)
     &   +(c4+c5*rr)*dexp(-c6*rr)*rr**2
c       move zero from asymptote to minimum
        v=v+4.74772265
        v=v/autoev

        dvdr=((c1+c2*rr)*(-c3)+c2)*dexp(-c3*rr)
     &      +((c4+c5*rr)*(-c6)+c5)*dexp(-c6*rr)*rr**2
     &       +(c4+c5*rr)*dexp(-c6*rr)*rr*2.d0
        dvdr=dvdr/autoev

        elseif (at(1).eq.29.and.at(2).eq.29) then
! O2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       Jasper April 3, 2012
        de=42046.5d0 ! exp De in cm-1 
        re=1.2075d0 ! exp in A
        c1= 2.6938139d0 ! my fit
        c2= 0.384763939d0
        c3= 0.812506485d0

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2    ! A and cm-1
        v=v/autocmi  ! convert to au

c        print *,rr,yy,beta,v

        dvdr=c1+2.d0*c2*yy+3.d0*c3*yy**2
        dvdr=dvdr*2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)
        dvdr=dvdr*autoang/autocmi  ! convert to au

c        print *,dvdr

        elseif (at(1).eq.26.and.at(2).eq.26) then
! N2 bath
!       fit to MRCI+Q/CBS(AQZ,A5Z) full valence
!       agrees reasonably well with more complicated form of LeRoy (JCP
!       125, 164310 (2006))
!       Jasper June 9, 2010
        de=79845.d0 ! exp De in cm-1 (Ronin, Luanay, Larzillier, 
!                                     PRL 53, 159 (1984), as quoted by
!                                     LeRoy)
        re=1.097679d0 ! exp in A
        c1=2.68872341 ! my fit
        c2=0.240070803
        c3=0.472261727

        yy=rr*autoang-re
        beta = c1+c2*yy+c3*yy**2
        v = de*(1.d0-dexp(-beta*yy))**2    ! A and cm-1
        v=v/autocmi  ! convert to au

c        print *,rr,yy,beta,v

        dvdr=c1+2.d0*c2*yy+3.d0*c3*yy**2
        dvdr=dvdr*2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)
        dvdr=dvdr*autoang/autocmi  ! convert to au

c        print *,dvdr

        elseif ((at(1).eq.27.and.at(2).eq.28).or.
     &          (at(1).eq.28.and.at(2).eq.27)) then
! CO bath
!       Morse. Fit to RKR data of PAUL H. KRUPENIE and STANLEY WEISSMAN,
!       J. CHEM. PHYS. 43, 1529 (1965)
!       with De = 11.06 eV
        de=11.06d0    ! exp De in eV
        de=de/autoev
        re=1.128322d0 ! exp in A
        re=re/autoang
        beta=1.d0/0.428d0  ! my fit in 1/A
        beta=beta*autoang

        yy=rr-re
        v = de*(1.d0-dexp(-beta*yy))**2
        dvdr=2.d0*de*(1.d0-dexp(-beta*yy))*dexp(-beta*yy)*beta

!       elseif (at(1).eq.??.and.at(2).eq.??) then
! OTHER DIATOMIC BATHS HERE
        else
          print *,"Don't know this diatomic bath"
          stop
        endif

      dvdx(1) =  dvdr*dx/rr
      dvdx(2) = -dvdr*dx/rr
      dvdy(1) =  dvdr*dy/rr
      dvdy(2) = -dvdr*dy/rr
      dvdz(1) =  dvdr*dz/rr
      dvdz(2) = -dvdr*dz/rr

      endif

      return
      end




c **********************************************************************
c **********************************************************************
      subroutine tinkerpot(symb,xx,yy,zz,pema,dxx,dyy,dzz,nat,mnat,idum)

      use sizes
      use atoms
      use files
      use inform
      use iounit

      implicit none

      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      double precision autokcal,autoang
      parameter(autokcal=627.509d0)
      parameter(autoang=0.52917706d0)

      integer iprepot,mnat,nat,i,idum,itinkxyz,itinkxyzlast
      double precision xx(mnat), yy(mnat), zz(mnat),
     &                dxx(mnat),dyy(mnat),dzz(mnat),
     &                pema,energy,derivs(3,maxatm)
      character*2 symb(mnat),dum

      integer nfitparams,mfitparams
      parameter (mfitparams=50)
      double precision ccc(mfitparams)
      common/amfit/ccc,nfitparams
      common/tinkerhack/itinkxyz

      save iprepot,itinkxyzlast

      entry prepot

      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

!      print *,"Got to prepot on proc ",my_id

      if (iprepot.eq.0) then
        call prepot2
        call prepot3
c     set up the structure and mechanics calculation
!        IF (my_id.eq.0) THEN
        call initial
        call getxyz
        call getkey
        call mechanic
!        ENDIF
        iprepot=1
        itinkxyz=0
        itinkxyzlast=0
!        call MPI_BCAST(nat, 1, MPI_INTEGER,
!     &               0, MPI_COMM_WORLD, ierr)
!        call MPI_BCAST(xx, n, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)
!        call MPI_BCAST(yy, n, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)
!        call MPI_BCAST(zz, n, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)
        return
      endif

      itinkxyz=idum
      if (itinkxyz.ne.itinkxyzlast) then
        call initial
        call getxyz
        call getkey
        call mechanic
      endif

      do i=1,n
        x(i)=xx(i)*autoang
        y(i)=yy(i)*autoang
        z(i)=zz(i)*autoang
      enddo

!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!      call MPI_BCAST(n, 1, MPI_INTEGER,
!     &               0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(x, n, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(y, n, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(z, n, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)

!      IF (my_id.eq.0) THEN
      call gradient (energy,derivs)
      pema=energy/autokcal

      do i=1,nat
        dxx(i)=derivs(1,i)/autokcal*autoang
        dyy(i)=derivs(2,i)/autokcal*autoang
        dzz(i)=derivs(3,i)/autokcal*autoang
      enddo
!      ENDIF

!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!      call MPI_BCAST(nat, 1, MPI_INTEGER,
!     &               0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(dxx, nat, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(dyy, nat, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)
!      call MPI_BCAST(dzz, nat, MPI_DOUBLE_PRECISION,
!     &               0, MPI_COMM_WORLD, ierr)

      itinkxyzlast=itinkxyz

      return

      end
c **********************************************************************
c **********************************************************************




********************************************

      subroutine prepot2

      implicit double precision(a-h,o-z)

      include 'mpif.h'
      integer my_id,nproc,ierr,term1,term2,ndata,t
      integer status(MPI_STATUS_SIZE)

      include 'paramls.inc'
      dimension iagroup(maxatom),
!     &  ind(maxterm,maxpair),
     &  nfound(maxatom),iatom(maxperm,maxatom),
     &  idum(maxatom),nngroup(maxatom),idum2(maxperm,maxatom),
     &  idum3(maxperm,maxatom),nperm0(maxperm),nperm1(maxperm),
     &  ibasis(maxterm),r(maxpair),
     &  basis(maxterm),dbasisdr(maxterm,maxpair),
     &  rrr(maxdata,maxpair),index(maxatom,maxatom),ix(maxperm,maxpair)
      double precision, allocatable :: basismpi(:), dbasisdrmpi(:,:)
!      double precision, allocatable :: basis(:),dbasisdr(:,:)
      integer, allocatable :: ind(:,:)
      character*2 symb(maxatom),dum
      logical lreadbasis
 
      common/foox/rrr,nncoef,natom1

      save npairs,nterms,ind,ibasis

      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

ccc GENERATE BASIS ccc
      if (.false.) then
ccc GENERATE BASIS ccc

c     generate atom permutation lists
      do i=1,natom
        nngroup(i)=0
      enddo
      do i=1,natom
        if (iagroup(i).gt.ngroup) ngroup=iagroup(i)
        nngroup(iagroup(i))=nngroup(iagroup(i))+1
      enddo

      nn=0

      do i=1,ngroup
        n=0
        do k=1,natom
          if (iagroup(k).eq.i) then
            n=n+1
            idum(n)=k
          endif
        enddo
      
        npermute=0
        call heapp(idum,n,n,idum2,npermute)
        nperm0(i)=nn+1
        nperm1(i)=nn+npermute
        do k=1,npermute
          nn=nn+1
          m=0
          do j=1,natom
            idum3(nn,j)=0
            if (iagroup(j).eq.i) then
              m=m+1
              idum3(nn,j)=idum2(k,m)
            endif
          enddo
        enddo
      enddo

      ntmp=1
      do i=1,ngroup
        idum(i)=nperm0(i)
        IF (my_id.eq.0)
     &    print *,
     &   "Group ",i," has ",(nperm1(i)-nperm0(i)+1)," permutations"
        ntmp=ntmp*(nperm1(i)-nperm0(i)+1)
      enddo
      IF (my_id.eq.0) print *,"For a total of ",ntmp," permutations"

      npermute=0
      do while (.true.)
        npermute=npermute+1
        if (npermute.gt.maxperm) then
          IF (my_id.eq.0) THEN
            print *,"npermute (",npermute,") > maxperm (",maxperm,")"
            print *,"NOTE: maxperm needs to be at least npermute + 1"
          ENDIF
          stop
        endif

        do i=1,natom
          iatom(npermute,i)=0
          do j=1,ngroup
            iatom(npermute,i)=iatom(npermute,i)+idum3(idum(j),i)
          enddo
        enddo

        idum(ngroup)=idum(ngroup)+1
 777    continue

        do i=1,ngroup
          if (idum(i).gt.nperm1(i)) then
            if (i.eq.1) go to 778
            idum(i)=nperm0(i)
            idum(i-1)=idum(i-1)+1
            go to 777
          endif
        enddo 

      enddo
 778  continue

      IF (my_id.eq.0) THEN
      write(6,*)
      write(6,*)'Atom permutations',npermute
      do i=1,min(npermute,100)
        write(6,*)i,":",(iatom(i,j),j=1,natom)
      enddo
      ENDIF

      ii=0
      do i=1,natom1
        do j=natom1+1,natom
          ii=ii+1
          index(i,j)=ii
        enddo
      enddo

      IF (my_id.eq.0) THEN
      write(6,*)
      write(6,*)"Pair permutations"
      write(6,'(22x,100(a3,"- ",a3,4x))')
     &   ((symb(i),symb(j),j=1+natom1,natom),i=1,natom1) 
      write(6,'(21x,100(i3," -",i3,4x))')((i,j,j=1+natom1,natom),
     &   i=1,natom1) 
      ENDIF

      do ii=1,npermute
        iix=0
        do i=1,natom1
          do j=natom1+1,natom
            iix=iix+1
            ix(ii,iix)=index(iatom(ii,i),iatom(ii,j))
          enddo
        enddo
        IF (my_id.eq.0) THEN
          if (ii.le.100) write(6,*)ii,":",(ix(ii,iix),iix=1,npairs)
        ENDIF
      enddo

c generate terms using individual power constraints
      allocate(ind(nterms,npairs))

      ii=1
      do i=1,npairs
        ind(ii,i)=0
      enddo
      do while (.true.)
        ii=ii+1
        if (ii.gt.maxterm) then
          IF (my_id.eq.0)
     &     print *,"number of terms (",ii,") > maxterm (",maxterm,")"
          stop
        endif

        do i=1,npairs
          ind(ii,i)=ind(ii-1,i)
        enddo
        ind(ii,npairs)=ind(ii,npairs)+1
 300    continue
        indtot=0
        do i=1,npairs
          indtot=indtot+ind(ii,i)
          ! ipow(i) would allow atom-atom-type-dependent limits
          if (ind(ii,i).gt.ipow.or.indtot.gt.ipowt) then
            if (i.eq.1) go to 400
            ind(ii,i)=0
            ind(ii,i-1)=ind(ii,i-1)+1
            go to 300
          endif
        enddo
      enddo
 400  continue
      nterms=ii-1

      IF (my_id.eq.0) THEN
       write(6,*) 
       write(6,*)"Basis # (Group):  Powers"
      ENDIF

c symmetrize
      nbasis=0
      DO ii=1,nterms
        ifail=0
        do i=1,ii-1
          do j=1,npermute
            ifail=1
            do k=1,npairs
              if (ind(i,k).ne.ind(ii,ix(j,k))) ifail=0
            enddo
            if (ifail.eq.1) go to 1010
          enddo
        enddo
 1010 continue

        if (ifail.eq.0) then
          nbasis=nbasis+1
          ibasis(ii)=nbasis
        else
          ibasis(ii)=ibasis(i)
        endif
        IF (my_id.eq.0)
     &   write(6,'(i5,"  (",i5,"):",100i8)')
     &     ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      ENDDO

      nncoef=nbasis
      IF (my_id.eq.0) THEN
      write(6,*)'nncoef = ',nncoef
 
      open(55,file="basis.dat")
      write(55,*)natom1,npairs,nncoef,nterms,
     & " ! atom pairs, coefficients, terms"
      write(55,*)" TERM GROUP :     EXPONENTS"
      do ii=1,nterms
        write(55,'(2i6," : ",1000i5)')
     &   ii,ibasis(ii),(ind(ii,j),j=1,npairs)
      enddo
      close(55)
      ENDIF

ccc READ BASIS ccc
      else
        IF (my_id.eq.0) THEN
        open(55,file="basis.dat")
        read(55,*)natom1,npairs,nncoef,nterms
        read(55,*)
        ENDIF

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(npairs, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nterms, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)

        allocate(ind(nterms,npairs))

        IF (my_id.eq.0) THEN
        do i=1,nterms
          read(55,*)k,ibasis(k),dum,(ind(k,j),j=1,npairs)
        enddo
        close(55)
        ENDIF

        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(natom1, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(nncoef, 1, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ibasis, nterms, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(ind, npairs*nterms, MPI_INTEGER,
     &                 0, MPI_COMM_WORLD, ierr)

!        print *,"npairs = ",npairs,"on proc",my_id
!        print *,"nterms = ",nterms,"on proc",my_id
!        print *,"ibasis(5) = ",ibasis(5),"on proc",my_id
!        print *,"ind(3,36) = ",ind(3,36),"on proc",my_id

      endif

      return

      entry funcs1(iii,basis,ncoef,dbasisdr)

      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(ncoef, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!      print *,"ncoef = ",ncoef,"in funcs1 on proc",my_id

      do i=1,ncoef
        basis(i)=0.d0
        do j=1,npairs
          dbasisdr(i,j)=0.d0
        enddo
      enddo

!      print *,"initialized basis and dbasisdr on proc",my_id

      do j=1,npairs
        r(j)=dexp(-rrr(iii,j)*autoang)
      enddo

      ndata=nterms/nproc+1*min(mod(nterms,nproc),1)

!      print *,"ndata = ",ndata,"on proc",my_id

      allocate(basismpi(nterms))
      allocate(dbasisdrmpi(nterms,npairs))

      do i=1,ncoef
        basismpi(i)=0.d0
        do j=1,npairs
          dbasisdrmpi(i,j)=0.d0
        enddo
      enddo
!      print *,"ind(1,j) =",(ind(1,j),j=1,npairs),"in funcs1"

ccccc MPI DO
      do i=my_id+1,nterms,nproc
        arg=1.d0
        do j=1,npairs
          if (ind(i,j).ne.0) arg=arg*(r(j)**ind(i,j))
        enddo
        basismpi(ibasis(i))=basismpi(ibasis(i))+arg
        do j=1,npairs
          dbasisdrmpi(ibasis(i),j)=dbasisdrmpi(ibasis(i),j)    ! dV/dy * dy/dr
     &      -arg*dble(ind(i,j))*autoang
        enddo
      enddo
ccccc MPI DO

      print *,"basis(1) =",basismpi(1),"in funcs1 on proc",my_id

      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

      call MPI_ALLREDUCE(basismpi, basis, nterms, MPI_DOUBLE_PRECISION,
     &                MPI_SUM, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE(dbasisdrmpi, dbasisdr, nterms*npairs,
     &                MPI_DOUBLE_PRECISION, MPI_SUM,
     &                MPI_COMM_WORLD, ierr)

      return 
      end

***************************************************

      recursive subroutine heapp(ia,size,n,iia,ii)

      include 'paramls.inc'
      integer i,n,size,ii
      integer ia(maxatom)
      integer iia(maxperm,maxatom)
      integer iagroup(maxatom)

      if (size.eq.1) then
         ii=ii+1
         do i=1,n
           iia(ii,i)=ia(i)
         enddo
        return
      endif

      do i=1,size
        call heapp(ia,size-1,n,iia,ii)
        if (mod(size,2).eq.1) then
          tmp=ia(1)
          ia(1)=ia(size)
          ia(size)=tmp
        else
          tmp=ia(i)
          ia(i)=ia(size)
          ia(size)=tmp
      endif

      enddo

      end subroutine

***************************************************
      subroutine lsfit(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatomx)
c
      implicit double precision (a-h,o-z)
c
      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)
c
      include 'paramls.inc'
      dimension coef(maxterm),sig(maxdata)
      dimension basis(maxterm)
      dimension dbasisdr(maxterm,maxpair),dvdr(maxpair)

      dimension vv(maxdata),rrr(maxdata,maxpair)
      dimension vv2(maxdata),rrr2(maxdata,maxpair)
      dimension rcom2(maxdata),xprint(50,20)
      dimension x(maxatom),y(maxatom),z(maxatom)

      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      integer at(maxatom),npairs

      logical lnumerical

      character*2 dum

      common/foox/rrr,nncoef,natom1

      save iprepot3
      entry prepot3

ccccc MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

      call MPI_BCAST(iprepot3, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
!      print *,"iprepot3 = ",iprepot3,"on proc",my_id

      if (iprepot3.eq.0) then
        IF (my_id.eq.0) THEN
        open(77,file="coef.dat")
        do k=1,nncoef
          read(77,*)i,coef(k)
        enddo
        ENDIF
        iprepot3=1
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        call MPI_BCAST(coef, nncoef, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
        return
      endif

c      print *
c      do j=1,natom
c        print *,at(j),x(j)*autoang,y(j)*autoang,z(j)*autoang
c      enddo

      ii=0
      do j=1,natom  ! molecule
      do k=j+1,natom  ! bath
        ii=ii+1
        rrr(1,ii)=dsqrt((x(j)-x(k))**2+(y(j)-y(k))**2+(z(j)-z(k))**2)
      enddo  
      enddo  

      ncoef=nncoef
      npairs=natom*(natom-1)/2

!      print *,"ncoef =",ncoef,"in prepot3 on proc",my_id
!      print *,"npairs =",npairs,"in prepot3 on proc",my_id
      
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call funcs1(1,basis,ncoef,dbasisdr) 
!      call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!      print *,"Got to 877 in prepot3 on proc",my_id
!      print *,"basis(1) =",basis(1),"in lsfit on proc",my_id

      IF (my_id.eq.0) THEN
      v=0.d0
      do j=1,ncoef
        v=v+coef(j)*basis(j)
      enddo
      v=v/autocmi
      ENDIF
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
      call MPI_BCAST(v, 1, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)

      do k=1,ii
        dvdr(k)=0.d0
        do j=1,ncoef
          dvdr(k)=dvdr(k)+coef(j)*dbasisdr(j,k)
        enddo
        dvdr(k)=dvdr(k)/autocmi
      enddo

      lnumerical=.true.
      lnumerical=.false.
      resp=0.0001
      ii=0
      do j=1,natom  ! molecule
      do k=j+1,natom  ! bath
        dx=x(j)-x(k)
        dy=y(j)-y(k)
        dz=z(j)-z(k)
        ii=ii+1

        if (lnumerical) then
          rrr(1,ii)=rrr(1,ii)+resp
          call MPI_BARRIER(MPI_COMM_WORLD, ier)
          call funcs1(1,basis,ncoef,dbasisdr) 
          vp=0.d0
          do l=1,ncoef
            vp=vp+coef(l)*basis(l)
          enddo
          vp=vp/autocmi
          rrr(1,ii)=rrr(1,ii)-2.d0*resp
          call MPI_BARRIER(MPI_COMM_WORLD, ier)
          call funcs1(1,basis,ncoef,dbasisdr) 
          vm=0.d0
          do l=1,ncoef
            vm=vm+coef(l)*basis(l)
          enddo
          vm=vm/autocmi
          rrr(1,ii)=rrr(1,ii)+resp
          dtmpdrr=(vp-vm)/(2.d0*resp)
        else
          dtmpdrr=dvdr(ii)
        endif
        dvdx(j) = dvdx(j) + dtmpdrr*dx/rrr(1,ii)
        dvdx(k) = dvdx(k) - dtmpdrr*dx/rrr(1,ii)
        dvdy(j) = dvdy(j) + dtmpdrr*dy/rrr(1,ii)
        dvdy(k) = dvdy(k) - dtmpdrr*dy/rrr(1,ii)
        dvdz(j) = dvdz(j) + dtmpdrr*dz/rrr(1,ii)
        dvdz(k) = dvdz(k) - dtmpdrr*dz/rrr(1,ii)
      enddo
      enddo

      return
 
      end

***************************************************

      subroutine rgexp(at,x,y,z,v,dvdx,dvdy,dvdz,natom,maxatom)

      implicit real*8(a-h,o-z)

      include 'mpif.h'
      integer my_id,nproc,ierr
      integer status(MPI_STATUS_SIZE)

      dimension x(maxatom),y(maxatom),z(maxatom)
      dimension dvdx(maxatom),dvdy(maxatom),dvdz(maxatom)
      parameter(autoev=27.2113961d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)
      parameter(autoang=0.529177249d0)
      integer at(maxatom)
      logical troya,cutoff

ccccc MPI
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_id, ierr)

      v1=0.d0
      v=0.d0
      do i=1,natom
        dvdz(i)=0.d0
        dvdy(i)=0.d0
        dvdx(i)=0.d0
      enddo

      do 1 i=1,natom
      do 2 j=i+1,natom
        m1=min(at(i),at(j))
        m2=max(at(i),at(j))
        troya=.false.   ! do or don't use Troya's form
        cutoff=.false.   ! do or don't use cutoff

        if (m1.ge.21) then ! two bath gases, skip this pair
          go to 2
        endif
        if (m2.le.20) then ! no bath gas, skip this pair
          go to 2
        endif

        if (m2.eq.21) then      ! He-
          if (m1.eq.1) then !    H
            aa=   3.5613416554821224
            bb=  0.29093368293887040
            cc=   1.5165871681294294
            rrc=   7.8982383462402739
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=   6.4297309502938740
            bb=  0.27714239680012726
            cc=   7.6654531794284928
            rrc=   5.6542253496969996
            cutoff=.true.
            aa=(10.d0**aa)
          elseif (m1.eq.3) then !    C
            aa=   7.1586261657152201
            bb=  0.22296437500485278
            cc=   6.0565445218329756
            rrc=   4.9366655126949777
            cutoff=.true.
            aa=(10.d0**aa)
          elseif (m1.eq.11) then !    C
            aa=   3.8786164447888081
            bb=  0.38199003107170471
            cc=   5.5398312746827889
            rrc=   5.8827497562148130
            cutoff=.true.
            aa=(10.d0**aa)
          elseif (m1.eq.12) then !    C
            aa=   7.7137538430550698
            bb=  0.33439442868091906
            cc=   1.8940308894700681
            rrc=   1.5424439934773555
            cutoff=.true.
            aa=(10.d0**aa)
          else
            write(6,*)"Cant find a Ne-? interaction",m1,m2
            stop
          endif
        elseif (m2.eq.23) then  ! Ar-
          if (m1.eq.1) then !    H
            aa=   6.7894331224729942
            bb=  0.25764780707510565
            cc=   6.3220103933812419
            rrc=   2.7862498144913403
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=   7.3353439234776987
            bb=  0.30185045895527779
            cc=   8.1275506825305328
            rrc=   1.5960859779572629
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a Ar-? interaction"
            stop
          endif
        elseif (m2.eq.24) then  ! Kr-
          if (m1.eq.1) then !    H
c       Troya's fit to CH4+Kr  (J. Phys. Chem. 110, 10834 (2006))
            aa=13754.02d0
            bb=3.238d0
            cc=-621.784d0
            troya=.true.
c         reoptimized for i-propyl + KR, QCISD(T)/CBS CC
c          aa=10.**4.19044465
c          bb=3.41540574
c          cc=-350.001526
c          troya=.true.
c         reoptimized for n-propyl + KR, QCISD(T)/CBS CC
c          aa=10.**4.27774285
c          bb=3.47558214
c          cc=-457.029939
c          troya=.true.
          elseif (m1.eq.2) then !    C
c       Troya's fit to CH4+Kr  (J. Phys. Chem. 110, 10834 (2006))
            aa=112927.4d0
            bb=3.520d0
            cc=-268.460d0
            troya=.true.
c         reoptimized for i-proply + KR, QCISD(T)/CBS CC
c           aa=10.**5.6093936
c           bb=3.91874142
c           cc=-1102.49641
c          troya=.true.
c         reoptimized for n-propyl + KR, QCISD(T)/CBS CC
c          aa=10.**5.10611591
c          bb=3.51407819
c          cc=-683.828242
c          troya=.true.
          else
            write(6,*)"Cant find a Kr-? interaction"
            stop
          endif
        elseif (m2.eq.25) then  ! H2-
          if (m1.eq.1) then !    H
c       aa=4.91209143  ! compromise radial, perpendicular fit
c       bb=0.386484268
c       cc=6.27030549
c       rrc=2.50001526
            aa=5.23924375 ! radial only fit
            bb=0.31757561
            cc=5.53859066
            rrc=2.62514115
c       aa=5.68126469 ! perp only
c       bb=0.277195654
c       cc=5.45313883
c       rrc=2.55272073
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
c       aa=5.1396527 ! compromise radial, perp fit
c       bb=0.428481094
c       cc=3.24997711
c       rrc=2.94143498
            aa=6.91601306 ! radial only fit
            bb=0.224359874
            cc=5.47520371
            rrc=2.46021912
c       aa=6.34962615 ! perp only
c       bb=0.289641102
c       cc=1.51567125
c       rrc=3.
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a H2-? interaction"
            stop
          endif
        elseif (m2.eq.29) then  ! O2-
          if (m1.eq.1) then !    H
            aa=6.04197821
            bb=0.288725394
            cc=5.99391156
            rrc=2.63286233
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=7.57335124
            bb=0.245461135
            cc=6.92338328
            rcc=2.03799554
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a N2-? interaction"
            stop
          endif
        elseif (m2.eq.26) then  ! N2-
          if (m1.eq.1) then !    H
            aa=   5.7800768359738344
            bb=  0.32909861021818337
            cc=   6.5168861989364117
            rrc=   2.3977632560101290
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=   7.7351792128122945
            bb=  0.24670372105706628
            cc=   5.7878800505168204
            rrc=   5.4064174131459506
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a N2-? interaction"
            stop
          endif
        elseif (m2.eq.27) then  ! C of a CO bath-
          if (m1.eq.1) then !    H
            aa=  6.67976928  ! radial O-in
            bb=  0.292242805
            cc=  0.738883633
            rrc=1.45118564
c         aa=7.12680441  ! compromise
c         bb=0.22606708
c         cc=6.56242561
c         rrc=3.28461562
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=  7.1884518 ! radial O-in
            bb=  0.238867763
            cc=  8.74794763
            rrc=2.50135807
c         aa=7.500412   ! compromise 
c         bb=0.276250801
c         cc=1.68080691
c         rrc=7.13440352
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a CO-? interaction"
            stop
          endif
        elseif (m2.eq.28) then  ! O of a CO bath-
          if (m1.eq.1) then !    H
            aa=  4.99996948 ! radial O-in
            bb=  0.413102512
            cc=  0.90664388
            rrc=8.74816126
c         aa=4.99996948  ! compromise
c         bb=0.398960234
c         cc=1.99975585
c         rrc=2.4582049
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=  6.9365215 ! radial O-in
            bb=  0.30321543
            cc=  8.98471023
            rrc=2.29291665
c         aa=5.99981689  ! compromise
c         bb=0.339911191
c         cc=8.74993133
c         rrc=3.00070193
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a CO-? interaction"
            stop
          endif

        elseif (m2.eq.30) then  ! O of a CO2 bath-
          if (m1.eq.1) then !    H
            aa=   5.7133727949153048     ! ch4+co2
            bb=  0.34664998206011610     
            cc=  0.99596395742398969     
            rrc=   3.5620687469056436 
c       aa=   5.6387658042445690     ! c2h6+co2
c       bb=  0.34522593901600213     
c       cc=   1.0625790534821908     
c       rrc=   5.6347689149186948 
c       aa=   5.9901516284627743     ! c3h8+co2
c       bb=  0.29905813248286556     
c       cc=  0.97208045779935071     
c       rrc=   5.9735371053325608     
c       aa=   5.9209847834754674         ! c4h10+co2
c       bb=  0.30248934638573333     
c       cc=  0.89604781961453317     
c       rrc=   5.0931276601671032  
c       aa=   6.0983532061220771     ! c6h14+co2
c       bb=  0.28379400791362602     
c       cc=   1.1563034568850479     
c       rrc=   5.6418613679327567 
c       aa=   6.2307699172033040      ! c3h7+co2
c       bb=  0.28078259439798375     
c       cc=  0.80652734130348946     
c       rrc=   5.6271488121243758  
c       aa=   5.8219126690642584     ! c3h6+co2
c       bb=  0.32890695356181021     
c       cc=   1.0272195531782671     
c       rrc=   5.2815980787991510     
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=   8.1053533306129744     ! ch4+co2
            bb=  0.22307556972273146     
            cc=   9.4218777816646462     
            rrc=   1.3779406000132564    
c       aa=   8.3648581603959684     ! c2h6+co2
c       bb=  0.20694020308422825     
c       cc=   8.5704724716637735     
c       rrc=   1.5343104837388377  
c       aa=   8.5979283801801607     ! c3h8+co2
c       bb=  0.18659348558531499     
c       cc=   8.0461824481527415     
c       rrc=   2.4289320686809606
c       aa=   9.2898529742370837         ! c4h10+co2
c       bb=  0.17399546261931620     
c       cc=   7.5561591615819070     
c       rrc=   1.2914664889326399 
c       aa=   8.3052207253534824     ! c6h14+co2
c       bb=  0.20020514580094537     
c       cc=   7.4568712008659173     
c       rrc=   1.3174723360629113 
c       aa=   8.7392765199184428      ! c3h7+co2
c       bb=  0.16779423263555787     
c       cc=   7.5032064289548970     
c       rrc=   2.2462085471905877  
c       aa=   8.4707004207312639     ! c3h6+co2
c       bb=  0.18922517562981259     
c       cc=   8.5051801750170721     
c       rrc=   2.3605094126181272     
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a CO-? interaction"
            stop
          endif

        elseif (m2.eq.31) then  ! C of a CO2 bath-
          if (m1.eq.1) then !    H
            aa=   6.8122495219351631     ! ch4+co2
            bb=  0.21132412081889679     
            cc=   2.2344589532813202     
            rrc=   10.029962896073105  
c       aa=   6.4741249460036743     ! c2h6+co2
c       bb=  0.24256987588076295     
c       cc=   2.2053547153059609     
c       rrc=   10.059309296212756  
c       aa=   6.3380281547942587     ! c3h8+co2
c       bb=  0.25676808207162849     
c       cc=   2.1762193294308858     
c       rrc=   13.841503485275682 
c       aa=   6.1308172927379401     ! c4h10+co2
c       bb=  0.27633983411544699     
c       cc=   2.0141563881384790     
c       rrc=   8.6865355116935365   
c       aa=   6.0713153595093292     ! c6h14+co2
c       bb=  0.28426927903312593     
c       cc=   1.7697915105668933     
c       rrc=   11.186847816404187  
c       aa=   6.3259484513222644      ! c3h7+co2
c       bb=  0.25357253979220440     
c       cc=   2.2358877046019270     
c       rrc=   14.194762521775969     
c       aa=   6.3293880234883453     ! c3h6+co2
c       bb=  0.24680269819476922     
c       cc=   2.2505154616178755     
c       rrc=   13.563189790079699 
            aa=(10.d0**aa)
            cutoff=.true.
          elseif (m1.eq.2) then !    C
            aa=   6.9509847863272913     ! ch4+co2
            bb=  0.30929037241689722     
            cc=   8.0012438281930862     
            rrc=   2.9929378971160046   
c       aa=   6.8966248282645832     ! c2h6+co2
c       bb=  0.30822559750532807     
c       cc=   8.1914097675564577     
c       rrc=   2.4305608130382788  
c       aa=   7.3101641244443316     ! c3h8+co2
c       bb=  0.27476119778451435     
c       cc=   8.0627887926346578     
c       rrc=   2.4072590749854239    
c       aa=   7.2409438615333048     ! c4h10+co2
c       bb=  0.27574561819503168     
c       cc=   8.2705165551441286     
c       rrc=   2.5895410360010476 
c       aa=   6.9786478824500158     ! c6h14+co2
c       bb=  0.29568615202742515     
c       cc=   8.3054763653913195     
c       rrc=   2.3518679436098902 
c       aa=   7.2963984404572129      ! c3h7+co2
c       bb=  0.27480489753048171     
c       cc=   8.0433127621772087     
c       rrc=   2.5001324683146540     
c       aa=   6.9860239890374540     ! c3h6+co2
c       bb=  0.29425624565678776     
c       cc=   7.8452859154569490     
c       rrc=   2.9125520650864503   
            aa=(10.d0**aa)
            cutoff=.true.
          else
            write(6,*)"Cant find a CO-? interaction"
            stop
          endif

        else
          write(6,*)"Cant find a ?-? interaction"
          stop
        endif

        dx=x(i)-x(j)
        dy=y(i)-y(j)
        dz=z(i)-z(j)
        rr=dsqrt(dx*dx+dy*dy+dz*dz)
        rra=rr*autoang

! NOTE CANNOT HAVE BOTH TROYA FORM AND CUTOFF FORM

      if (troya) then   ! Troya uses different form & units
        v=aa*dexp(-rra*bb)+cc/rra**6
        v=v/autokcal
        dvdr = -aa*bb*dexp(-rra*bb)-6.d0*cc/rra**7
        dvdr=dvdr/autokcal*autoang
      elseif (cutoff) then  ! cutoff 1/R**-6 at short distances
        v=aa*dexp(-rra/bb)-(cc**6/(rra**6+rrc**6))
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)
     &      +6.d0*(cc**6)*(rra**5)/(rra**6+rrc**6)**2
        dvdr=dvdr/autocmi*autoang
      else
        v=aa*dexp(-rra/bb)-(cc/rra)**6
        v=v/autocmi
        dvdr = -aa/bb*dexp(-rra/bb)+(6.d0/rra)*(cc/rra)**6
        dvdr=dvdr/autocmi*autoang
      endif
      v1=v1+v

c      print *,m1,m2,rra,v*autocmi,v1*autocmi

c derivs = sum over all bonds (DV/DRij * DRij/DXi = DV/DRij * (Xi-Xj)/Rij)
      dvdx(i) = dvdx(i) + dvdr*dx/rr
      dvdx(j) = dvdx(j) - dvdr*dx/rr
      dvdy(i) = dvdy(i) + dvdr*dy/rr
      dvdy(j) = dvdy(j) - dvdr*dy/rr
      dvdz(i) = dvdz(i) + dvdr*dz/rr
      dvdz(j) = dvdz(j) - dvdr*dz/rr

    2 continue
    1 continue
      v=v1
c      print *,'rgexp',v1

      return

      end

c **********************************************************************

