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


      subroutine checkhop(phop,nsurf,lhop,newsurf)

c Used for surface hopping methods.
c Draw a random number and compare it to PHOP to check for a hop.
c If a hop occurs, LHOP = .TRUE. and NEWSURF = whichever surface
c the system is trying to hop to.  The hop does not actually 
c occur in this subroutine.  Hops (frustrated and successful) 
c occur in HOP.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

c input
      double precision phop(mnsurf)
      integer nsurf

c output
      logical lhop
      integer newsurf

c local
      integer i
      double precision dxi,tmp

c get random number
      dxi = sprng(rng_stream)

c check for hop
      newsurf = 0
      tmp = 0.d0
      do i=1,nsurft
      tmp = tmp + phop(i)
      if (dxi.le.tmp.and.newsurf.eq.0) newsurf = i
      enddo
      if (newsurf.eq.0) then
        write(6,*)"NEWSURF = 0 in HOPCHECK!"
        stop
      endif

      if (newsurf.ne.nsurf) then
         lhop = .true.
      endif

      return

      end
