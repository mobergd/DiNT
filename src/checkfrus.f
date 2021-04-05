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


      subroutine checkfrus(nclu,nsurf,pp,pem,dvec,frusflag,tu_maxt)

c Used for the FSTU method.
c At every step, check to see of a hop would be energetically allowed
c to all other states.  If a hop is allowed to state K, set FRUSFLAG(K) = 0.
c Otherwise, FRUSFLAG(K) = 1.

      implicit none
      include 'param.f'
      include 'c_sys.f'

      integer nclu,frusflag(mnsurf)
      integer nsurf

      integer is,i,j
      double precision pem(mnsurf),gpem(3,mnat,mnsurf),
     & dvec(3,mnat,mnsurf,mnsurf),pp(3,mnat),
     & sectrm,a,b,a2,ab,radi,scr,deltae,tu_maxt(mnsurf),
     & talongd

c check for forbidden hops
      do is = 1,nsurft
       if (is.ne.nsurf) then
        deltae = pem(is) - pem(nsurf)
        a=0.0d0
        b=0.0d0
        do i=1,3
        do j=1,nclu
         scr=dvec(i,j,nsurf,is)/mm(j)
c        a = momentum.dot.dvec/mass
c        b = dvec.dot.dvec/mass
         a=a+scr*pp(i,j)
         b=b+scr*dvec(i,j,nsurf,is)
        enddo
        enddo
        a2=a*a
        ab=a/b
c       sectrm = deltae/kinetic energy along d
        if (b.gt.0.d0) then
c       DVEC is nonzero
        sectrm = 2.0d0*deltae*b/a2
        talongd = a2/b/2.d0 ! kinetic energy along d
        tu_maxt(is) = 0.5d0/(deltae - talongd)
        radi = 1.0d0-sectrm
        else
c       DVEC is zero, mark as frustrated
        radi = -1.d0
        tu_maxt(is) = 0.d0
        endif
c       if radi > 0 then delta > kinetic energy along d, so we can hop
c       note that for a hop down, sectrm is negative, so we can always hop
        if(radi.ge.0.0d0) then
c        successful hop
         frusflag(is) = 0
        else
c        frustrated hop
         frusflag(is) = 1
        endif
       else
        frusflag(is) = 0
       endif
      enddo

      return
      end
