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


      subroutine ranclu(xx,nclu)

c Generate a random cluster of atoms within radius RADIUS,
c and with minimum and maximum minimum-atom-atom distances
c of RSMALL and RBIG, respectively.

      implicit none

      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer nclu,iclu,ikeep,itry,i,k,irej
      double precision xx(3,mnat),xyz(3),rr,rad
      double precision radius,rsmall,rbig,rrm
      real ran(3)

      radius = rran(1)
      rsmall = rran(2)
      rbig = rran(3)
      iclu = 0

 10   continue

      ikeep = 0
      iclu = iclu + 1
      if (iclu.gt.20) then
       write(6,*)"Tried too many times... quitting"
       stop
      endif
c      write(6,*)"Attempt #",iclu

c generate nclu points within the sphere
      itry = 0
      do while (ikeep < nclu)
        itry = itry + 1
        if (itry.eq.1000000) then
          write(6,*)"Bad packing.  Start again."
          go to 10
        endif
        do i=1,3
        ran(i) = sprng(rng_stream)
        enddo
        rad=0.d0
        do i=1,3
          xyz(i) = (2.d0*dble(ran(i))-1.d0)*radius
          rad = rad + xyz(i)**2
        enddo
        rad=dsqrt(rad)
        if (rad.le.radius) then
          irej = 0
          rrm = radius*100.d0
          if (ikeep.ne.0) then
            do i=1,ikeep
              rr = 0.d0
              do k=1,3
                rr = rr + (xyz(k) - xx(k,i))**2
              enddo
              rr=dsqrt(rr)
              rrm=min(rr,rrm)
            enddo
            if (rrm.lt.rsmall) irej = 1
            if (rrm.gt.rbig) irej = 1
          endif
          if (irej.eq.0) then
            ikeep = ikeep + 1
            do i=1,3
              xx(i,ikeep) = xyz(i)
            enddo
c            write(6,*)"keeping ",(xx(k,ikeep),k=1,3)
c            write(6,*)"keeping ",(xx(k,ikeep)*autoang,k=1,3)
          endif
        endif
      enddo

      return

      end

