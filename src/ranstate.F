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


      subroutine ranstate(freq,nmqn,temp)

      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>

c local
      integer nmax,i
      double precision temp,freq,ptot,ei,pp,p
      double precision r1,nmqn

      nmax=100
      ptot=0.d0
      do i=0,nmax
        ei=freq*(dble(i)+0.5d0)
        ptot=ptot+dexp(-ei/(temp*kb))
      enddo

      r1 = sprng(rng_stream)

      pp=0.d0
      nmqn=-1.d0
      do i=0,nmax
        ei=freq*(dble(i)+0.5d0)
        pp=pp+dexp(-ei/(temp*kb))
        p = pp/ptot
        if (r1.lt.p.and.nmqn.lt.0.d0) nmqn = dble(i)
      enddo
      if (nmqn.lt.0.d0) then
c     something went wrong. this can happen when the temperature is small and all the boltzmann weights are ~zero.
      nmqn=0.d0
      endif

      return

      end
