c
c   Dint – version 2.0  is licensed under the Apache License, Version 2.0 (the "License");
c   you may not use Dint – version 2.0 except in compliance with the License.
c   You may obtain a copy of the License at
c       http://www.apache.org/licenses/LICENSE-2.0
c   The license is also given in the LICENSE file.
c   Unless required by applicable law or agreed to in writing, software
c   distributed under the License is distributed on an "AS IS" BASIS,
c   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
c   See the License for the specific language governing permissions and limitations under the License.
c
c -------------------------------------------------------------------------------------------
c  Dint : Direct Nonadiabatic Trajectories A code for non-Born–Oppenheimer molecular dynamics 
c  
c  version 2.0                                    
c
c  A. W. Jasper                  
c  Argonne National Laboratories     
c
c  Rui Ming Zhang                 
c  Tsinghua University
c               
c  and                  
c    
c  D. G. Truhlar                 
c  University of Minnesota
c
c  copyright  2020
c  Donald G. Truhalar and Regents of the University of Minnesota 
c----------------------------------------------------------------------------------------------


      subroutine nmturn(symb,xx0,mm,nclu,vec,etot,rin,rout,
     & repflag,nsurf,nsurft)

c A VERSION OF RTBIS FROM NR

      implicit none
      include 'param.f'
      integer i,j,k,nclu,repflag,nsurf,nn,itype(3*mnat),mcpar(2,3*mnat),
     & nsurft
      double precision xx(3,mnat),vec(3,mnat),etot,rin,rout,
     & xx0(3,mnat),dd,mm(mnat),pema(mnsurf),pemd(mnsurf,mnsurf),
     & gpema(3,mnat,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),v(nsurf),x1,x2,f,fmid,root,dx,
     & xmid,dxtol
      character*2 symb(mnat)
      logical eprint

c     zero these for now
      do i=1,3*mnat
      itype(i)=0
      mcpar(1,i)=0
      mcpar(2,i)=0
      enddo

c     print energies while searching? 
      eprint=.true.
c      eprint=.false.

c     converge rin and rout to this tolerance
      dxtol = 1.d-3

c     outer turning point first, then inner below
      x1=1.  ! guesses for outer
      x2=0.

      nn=1

 12   call nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,
     & x1,itype,mcpar,v)
      f=v(nsurf)-etot

      call nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,
     & x2,itype,mcpar,v)
      fmid=v(nsurf)-etot

      if (eprint) write(6,*)f*autoev,fmid*autoev,x1,x2

      if (f*fmid.ge.0.d0) then
c        write(6,*)"root not bracketed in TURNNM"
        x1=x1*1.5
        go to 12
      endif
      if (f.lt.0.d0) then
          root=x1
          dx=x2-x1
      else
          root=x2
          dx=x1-x2
      end if
      do 11 j = 1, 500
          dx=dx*0.5d0
          xmid=root+dx

          call nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,
     & xmid,itype,mcpar,v)
          fmid= v(nsurf)-etot
          if (eprint) write(6,*)fmid*autoev,xmid

          if (fmid.le.0.d0) root=xmid
          if (abs(dx).lt.dxtol.or.fmid.eq.0.d0) then
            rout=root
            go to 20
          endif
   11 continue
      write(6,*)"Too many iterations in TURNNM"
      stop

c     inner turning point
 20   x1=-1.5 ! guesses for inner
      x2=0.

c      write(6,*)'found it'
 22   call nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,
     & x1,itype,mcpar,v)
      f=v(nsurf)-etot

      call nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,
     & x2,itype,mcpar,v)
      fmid=v(nsurf)-etot

      if (eprint) write(6,*)f*autoev,fmid*autoev,x1,x2

      if (f*fmid.ge.0.d0) then
c        write(6,*)"root not bracketed in TURNNM"
        x1=x1*1.5
        go to 22
        stop
      endif
      if (f.lt.0.d0) then
          root=x1
          dx=x2-x1
      else
          root=x2
          dx=x1-x2
      end if
      do 21 j = 1, 500
          dx=dx*0.5d0
          xmid=root+dx

          call nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,xmid,
     & itype,mcpar,v)
          fmid= v(nsurf)-etot
          if (eprint) write(6,*)fmid*autoev,xmid

          if (fmid.le.0.d0) root=xmid
          if (abs(dx).lt.dxtol.or.fmid.eq.0.d0) then
            rin=root
            return
          endif
   21 continue
      write(6,*)"Too many iterations in TURNNM"
      stop










      end

