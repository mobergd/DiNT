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

      subroutine popnorm(xxm,ppm,mmm,nclu,im)

c Populate the normal modes using information precalculated in
c NORMOD.  This subroutine selects initial conditions when 
c INITx = 2.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
#include <sprng_f.h>

      integer i,j,nclu,ido,im,nbound
      double precision etmp,mmm(mnat),xxm(3,mnat),ppm(3,mnat),
     & xvec(mnat*3),sgn(mnat*3),xxnm(3,mnat),vvnm(3,mnat),
     & kinho,potho,displacement,rturntmp,tmp

c generate two sets of 3*N-6 random numbers
      if (nmtype(im).eq.0) then
        nbound=3*nclu-6
      else
        nbound=3*nclu-7
      endif
      if (nclu.eq.2) then
        nbound=1
      endif
      do i=1,nbound
        xvec(i) = sprng(rng_stream)
        sgn(i) = sprng(rng_stream)
        if (sgn(i).ge.0.5d0) then
          sgn(i) = 1.d0
        else
          sgn(i) = -1.d0
        endif
      enddo

c initialize displacement and velocity
      do i=1,3
      do j=1,nclu
        xxnm(i,j) = 0.d0
        vvnm(i,j) = 0.d0
      enddo
      enddo

c populate all modes as indicated
c NOTE: xxnm, rturntmp are mass-scaled
c freq is independent of mass-scaling
      do i=1,3
      do j=1,nclu
        ppm(i,j) = 0.d0
        do ido=1,nbound
        etmp = freq(ido,im)*(0.5d0+nmqn(ido,im))
        rturntmp = dsqrt(2.d0*etmp/(freq(ido,im)**2*mu))
        displacement = rturntmp*dcos(pi*xvec(ido))
c       project displacement onto normal mode vector (mass-scaled)
c       sum displacements for all modes
        xxnm(i,j) = xxnm(i,j) + nmvec(i,j,ido,im)*displacement
c          write(6,*)i,j,ido,im,nmvec(i,j,ido,im)
c       compute HO estimate of potential
        potho = 0.5d0*freq(ido,im)**2*mu*displacement**2
c       the available kinetic energy (in the HO approx) is therefore
        kinho = max(0.d0,(etmp-potho))
c       assign velocity based on this energy
        vvnm(i,j) = vvnm(i,j)
     &    + sgn(ido)*nmvec(i,j,ido,im)*dsqrt(2.d0*kinho/mu)
        enddo
      enddo
      enddo

      if (nmtype(im).eq.1) then
      do i=1,3
      do j=1,nclu
        kinho = dabs(nmqn(nbound+1,im)) 
        tmp=1.d0
        if (nmqn(nbound+1,im).lt.0.d0) tmp=-1.d0
        vvnm(i,j) = vvnm(i,j)
     &    + tmp*nmvec(i,j,nbound+1,im)*dsqrt(2.d0*kinho/mu)
      enddo
      enddo
      endif

c un-scale and transform velocity to momentum
      do i=1,3
      do j=1,nclu
        xxm(i,j) = xxm(i,j) + xxnm(i,j)*dsqrt(mu/mmm(j))
        ppm(i,j) = ppm(i,j) + vvnm(i,j)*dsqrt(mu/mmm(j))*mmm(j)
c       write(6,*)i,j,ppm(i,j)
      enddo
      enddo


      return
      end
