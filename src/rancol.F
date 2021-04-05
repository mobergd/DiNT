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


      subroutine rancol(temp0,ecol)

c Generate collision energies from a Maxwell distribution.
c The probability distribution is
c    P(E)dE = E/kT exp(-E/kT) * dE/kT

c This is taken from gamdev from Numerical Recipes
c We select from 
c    p_a(x)dx = x exp(-x) dx/gamma(2)  (i.e., ia = 2)
c and set x = E/kT to obtain E

      implicit none
      include 'param.f'
c      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

c in/out
      double precision temp0,ecol

c local
      double precision x,r
      integer ia,i

      ia = 2
      x=1.d0
      do i=1,ia
         r=sprng(rng_stream)
         x=x*r
      enddo
      x=-dlog(x)

      ecol=x*kb*temp0   ! in au

      return

      end
