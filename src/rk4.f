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

      subroutine rk4(y,dydx,nv,hstep,nsurf)
c     Fourth order Runge-Kutta integrator. From numerical recipes.

      implicit none
      include 'param.f'

      integer nv
      double precision y(mnyarray),yt(mnyarray),dydx(mnyarray)

      integer i,j,nsurf
      double precision hh,h6,dum,hstep
      double precision dym(mnyarray),dyt(mnyarray)

      hh=hstep/2.d0
      h6=hstep/6.d0
      dum=0.d0

c     first step
      do 11 i=1,nv
      yt(i)=y(i)+hh*dydx(i)
 11   continue

c     second step
      call derivs(dum,yt,dyt,nv,nsurf)
      do 12 i=1,nv
      yt(i)=y(i)+hh*dyt(i)
 12   continue

c     third step
      call derivs(dum,yt,dym,nv,nsurf)
      do 13 i=1,nv
      yt(i)=y(i)+hstep*dym(i)
      dym(i)=dyt(i)+dym(i)
 13   continue

c     fourth step
      call derivs(dum,yt,dyt,nv,nsurf)
c     accumulate increments with proper weights
      do 14 i=1,nv
      y(i)=y(i)+h6*(dydx(i)+dyt(i)+2.d0*dym(i))
 14   continue

      return
      end
