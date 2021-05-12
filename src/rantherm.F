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
c  Donald G. Truhlar and Regents of the University of Minnesota 
c----------------------------------------------------------------------------------------------


      subroutine rantherm(pp,mmm,nclu,temp0)

c Generate momenta randomly according to a Maxwell distribution.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer nclu
      double precision pp(3,mnat),mmm(mnat),temp0

c local
      integer i,ix,iy,iz,j
      real ran(2)
      double precision r1,r2,x1,x2,w,y1,y2

c convert uniformly distributed numbers to a gaussian distributed set
c polar form of the Box-Muller transformation
c transformed random numbers are from a distribution that has zero mean
c and unit standard deviation
      do i=1,nclu
      do j=1,3
 10      continue
         r1 = sprng(rng_stream)
         r2 = sprng(rng_stream)
         x1=2.d0*r1-1.d0
         x2=2.d0*r2-1.d0
         w = x1*x1+x2*x2
         if (w.ge.1.d0.or.w.eq.0.d0) go to 10
         w = dsqrt(-2.d0*dlog(w)/w)
         y1 = x1*w
         y2 = x2*w
c you automatically get two random numbers (y1 and y2)
c i'm just going to use one for now, inefficient but simpler
c        transform to velocity kb*tempt, mmm, and pp are in a.u.
         pp(j,i)=y1*dsqrt( kb * temp0 * mmm(i) )
      enddo
      enddo

      return

      end
