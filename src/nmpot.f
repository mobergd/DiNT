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
c  Argonne National Laboratory     
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


      subroutine nmpot(symb,xx0,mm,nclu,nn,vec,repflag,nsurft,d,
     & itype,mcpar,v)

c     returns the energy EE of a structure displaced from XX0 along NN
c     mass-scaled vector VEC by some distance DD

      implicit none
      include 'param.f'
      integer i,j,k,nclu,repflag,nn,itype(3*mnat),mcpar(2,3*mnat),
     & a1,a2,n1,nsurft,maxstep,im
      double precision xx(3,mnat),vec(3,mnat,3*mnat),
     & xx0(3,mnat),d(3*mnat),mm(mnat),pema(mnsurf),pemd(mnsurf,mnsurf),
     & gpema(3,mnat,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),v(mnsurf),phi,gtmp,gtmp1,gapgrad,h12
      character*2 symb(mnat)

      common/lz/gapgrad,h12

      do j=1,3
      do k=1,nclu
        xx(j,k)=xx0(j,k)
      enddo
      enddo

c Do Cartesian normal mode displacements first
      do i=1,nn
        if (itype(i).eq.0) then ! Cartesian
        do j=1,3
        do k=1,nclu
          xx(j,k)=xx(j,k)+vec(j,k,i)*d(i)*dsqrt(mu/mm(k))
        enddo
        enddo
        endif
      enddo

c           write(80,*)nclu
c          write(80,1081)i,phi
c          do k=1,nclu
c            write(80,1080)symb(k),(xx(j,k)*autoang,j=1,3)
c          enddo

C Do torsions next
      do i=1,nn
        if (itype(i).eq.1) then ! torsion
        a1=mcpar(1,i)
        a2=mcpar(2,i)
        n1=a2-1
        phi=d(i)
           call torphi(xx,nclu,a1,a2,phi,n1)
c           write(80,*)nclu
c           write(80,1081)i,phi
c           do k=1,nclu
c             write(80,1080)symb(k),(xx(j,k)*autoang,j=1,3)
c           enddo
c 1080       format(a2,3f15.5)
c 1081       format(i7,50f15.5)
c          stop
        endif
      enddo

c get energy
      call getpem(xx,nclu,pema,pemd,gpema,gpemd,dvec,symb)
      do i=1,nsurft
      if (repflag.eq.1) then
        v(i) = pemd(i,i)
      else
        v(i) = pema(i)
      endif
      enddo

c          write(80,*)nclu
c          write(80,1081)0,v(1)*autoev,v(2)*autoev,gtmp,
c     & pemd(1,2)*autocmi
c          do i=1,nclu
c            write(80,1080)symb(i),(xx(j,i)*autoang,j=1,3)
c 1080       format(a2,3f15.5)
c 1081       format(i7,50f15.5)
c          enddo

c      gapgrad=gtmp
c      h12=pemd(1,2)

      return
      end

