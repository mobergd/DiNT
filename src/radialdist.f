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


      subroutine radialdist(itmp,xx,nclu,step,time,
     &   nbinrad,raddist,iprint)

c really the pair correlation function

      implicit none
      include 'param.f'

      integer i,j,k,nclu,ibin,nbinrad,iprint,itmp
      double precision xx(3,mnat),rij,raddist(0:nbinrad+1),tmp,
     & step,rbin,time,xnorm,rmax,rmin

      rmin = 0.d0 / autoang
      rmax = 30.d0 / autoang
      rbin = (rmax-rmin)/dble(nbinrad)
    
      do i=1,nclu
      do j=i+1,nclu
        rij = 0.d0
        do k=1,3
        rij = rij + (xx(k,i)-xx(k,j))**2
        enddo
        rij = dsqrt(rij)
        if (rij.le.rmin) then
          ibin = 0
        elseif (rij.ge.rmax) then
          ibin = nbinrad + 1
        else
          rij = rij - rmin
          ibin = int(rij/rbin) + 1 ! determine bin number
          ibin = max(1,ibin)
          ibin = min(nbinrad,ibin)
        endif
        raddist(ibin) = raddist(ibin) + step
      enddo
      enddo

      if (iprint.eq.1) then
      tmp = time*dble(nclu)
      tmp = tmp*rbin*4.d0*pi
c      http://www.physics.emory.edu/~weeks/idl/gofr2.html
      write(41,10)itmp,time*autofs,raddist(0)/time,
     &  (raddist(k)/(tmp*((dble(k)-0.5d0)*rbin+rmin)**2),k=1,nbinrad),
     &   raddist(nbinrad+1)/time
      endif
 10   format(i7,100f12.4)

      return
 
      end
