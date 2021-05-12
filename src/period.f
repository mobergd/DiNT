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


      subroutine period(edia,jj,arr,mm,rin,rout,peri,nsurf)

c used by the atom-diatom initial conditions method
c computes the period in au of the diatomic vibration.

      implicit none
      include 'param.f'

c input
      integer arr,nsurf,jj
      double precision rin,rout,mm(mnat),edia,xj

c output
      double precision peri

c local
      integer ncheb,i
      double precision arg,x(12),rmid,rdif,r,term,v,sum,rmass

      xj = dble(jj)

      ncheb = 12
      do 5 i = 1,ncheb
        arg = dble(2*i-1)/dble(2*ncheb)*pi
        x(i) = dcos(arg)
    5 continue

      rdif = (rout-rin)/2.d0
      rmid = (rout+rin)/2.d0
      sum = 0.d0
      do 10 i = 1,ncheb
        r = rmid+rdif*x(i)
        call diapot(r,arr,v,xj,mm,nsurf)
        term = sqrt((r-rin)*(rout-r)/(edia-v))
        sum = sum+term
   10 continue
      if (arr.eq.1) rmass = mm(1)*mm(2)/(mm(1)+mm(2))
      if (arr.eq.2) rmass = mm(2)*mm(3)/(mm(2)+mm(3))
      if (arr.eq.3) rmass = mm(3)*mm(1)/(mm(3)+mm(1))
      peri = sqrt(2.d0*rmass)*sum*pi/dble(ncheb)

      return

      end

