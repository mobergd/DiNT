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


      subroutine ransamp(lbins,ifems,ntot,filex2,filep2,nn,xx,pp,molpe)


      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>

c in/out
      integer ntot,nn
      character*10 filex2,filep2
      character*16 filexx,filepp
      character*2 dum
      double precision pp(3,mnat),xx(3,mnat),molpe
      logical lbins,ifems

c local
      integer n,i,j,nrcl,itmp1,itmp2, istep
      double precision r
      character*2 ctmp(nn)

      r=sprng(rng_stream)
      r=r*dble(ntot)
      n=nint(r)

c      filexx="../r0/" // filex2
c      filepp="../r0/" // filep2
      filexx=filex2
      filepp=filep2

      write(6,*)"Reading sample ",n

      if (lbins) then
 1080   format(a2,3f15.5)
 1081   format(i5,i7,f15.5)
        nrcl=(3+5*nn)*10
        open(93,file=filexx,form='unformatted',
     &   status='unknown',access='direct',RECL=nrcl)
        read(93,REC=n)itmp1,itmp2,molpe,(ctmp(i),i=1,nn),
     &   ((xx(j,i),j=1,3),i=1,nn)
        write(6,1081)itmp1,itmp2,molpe
        do i=1,nn
          write(6,1080)ctmp(i),(xx(j,i),j=1,3)
          xx(1,i) = xx(1,i)/autoang
          xx(2,i) = xx(2,i)/autoang
          xx(3,i) = xx(3,i)/autoang
        enddo
      else
        open(unit=93,file=filexx)
        read(93,*)
        read(93,*)
        do i=1,n-1
          read(93,*)
          read(93,*)
          do j=1,nn
            read(93,*)
          enddo
        enddo
        read(93,*)
        read(93,*)istep,molpe
        do j=1,nn
          read(93,*)dum,xx(1,j),xx(2,j),xx(3,j)
          xx(1,j)=xx(1,j)/autoang
          xx(2,j)=xx(2,j)/autoang
          xx(3,j)=xx(3,j)/autoang
        enddo
      endif
      close(93)
      
      if (ifems) return
      if (lbins) then
        open(94,file=filepp,form='unformatted',
     &   status='unknown',access='direct',RECL=nrcl)
        read(94,REC=n)itmp1,itmp2,molpe,(ctmp(i),i=1,nn),
     &   ((pp(j,i),j=1,3),i=1,nn)
        write(6,1081)itmp1,itmp2,molpe
        do i=1,nn
          write(6,1080)ctmp(i),(pp(j,i),j=1,3)
        enddo
      else
        open(unit=94,file=filepp)
        read(94,*)
        read(94,*)
        do i=1,n-1
          read(94,*)
          read(94,*)
          do j=1,nn
            read(94,*)
          enddo
        enddo
        read(94,*)
        read(94,*)istep,molpe
        do j=1,nn
          read(94,*)dum,pp(1,j),pp(2,j),pp(3,j)
        enddo
      endif
      close(94)

      return

      end
