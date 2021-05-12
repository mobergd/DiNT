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


      subroutine xptoy(xx,pp,cre,cim,gv,gcre,gcim,mm,y,dydx,
     &   nat,nsurft,ido,methflag,pem,phase)

c The integrators integrate a single array Y using its derivatives
c DYDX.  This subroutine transforms back and forth between the 
c Y array and the various quantities.
c The Y array contains:
c   Elements 1 to NAT           -> X coordinates of all atoms XX(1,i)
c   Elements NAT+1 to 2*NAT     -> Y coordinates of all atoms XX(2,i)
c   Elements 2*NAT+1 to 3*NAT   -> Z coordinates of all atoms XX(3,i)
c   Elements 3*NAT+1 to 4*NAT   -> X momenta of all atoms PP(1,i)
c   Elements 4*NAT+1 to 5*NAT   -> Y momenta of all atoms PP(2,i)
c   Elements 5*NAT+1 to 6*NAT   -> Z momenta of all atoms PP(3,i)
c   Elements 6*NAT+1 to 6*NAT+NSURFT   -> PHASE(k)
c   Elements 6*NAT+NSURFT+1 to 6*NAT+2*NSURFT   -> Real part of the electronic
c                                                  variables CRE(k)
c   Elements 6*NAT+2*NSURFT+1 to 6*NAT+3*NSURFT -> Imaginary part of the electronic
c                                                  variables CRE(k)
c   Elements 6*NAT+3*NSURFT+1 to 6*NAT+4*NSURFT -> Real part of the electronic
c                                                  variables (coherent part for CSDM)
c   Elements 6*NAT+4*NSURFT+1 to 6*NAT+5*NSURFT -> Imaginary part of the electronic
c                                                  variables (coherent part for CSDM)

      implicit none
      include 'param.f'

      integer i,j,ido,i2,methflag
      integer nat,nsurft
      double precision xx(3,mnat),pp(3,mnat),gv(3,mnat),mm(mnat)
      double precision y(mnyarray),dydx(mnyarray),
     &  cre(2*mnsurf),cim(2*mnsurf),gcre(2*mnsurf),gcim(2*mnsurf),
     &  pem(mnsurf),phase(mnsurf)

      if (ido.eq.1) then

c     break y into components
c     get info from y and dydx
      do i=1,nat
      do j=1,3
      xx(j,i)=y((j-1)*nat+i)
      pp(j,i)=y((j-1)*nat+i+3*nat)
      gv(j,i)=-dydx((j-1)*nat+i+3*nat)
      enddo
      enddo
      do i=1,nsurft
      phase(i) = y(6*nat+i)
      enddo
      do i=1,nsurft
      cre(i) = y(6*nat+nsurft+i)
      cim(i) = y(6*nat+2*nsurft+i)
      gcre(i) = dydx(6*nat+nsurft+i)
      gcim(i) = dydx(6*nat+2*nsurft+i)
      enddo

      if (methflag.eq.4) then
c     CSDM (coherent terms)
      do i=1,nsurft
      i2 = i + nsurft
      cre(i2) = y(6*nat+3*nsurft+i)
      cim(i2) = y(6*nat+4*nsurft+i)
      gcre(i2) = dydx(6*nat+3*nsurft+i)
      gcim(i2) = dydx(6*nat+4*nsurft+i)
      enddo
      endif

      else
c     assemble components into a 1-D vector for integration
      do i=1,nat
      do j=1,3
      y((j-1)*nat+i)=xx(j,i)
      y((j-1)*nat+i+3*nat)=pp(j,i)
      dydx((j-1)*nat+i)=pp(j,i)/mm(i)
      dydx((j-1)*nat+i+3*nat)=-gv(j,i)
      enddo
      enddo
      do i=1,nsurft
      y(6*nat+i) = phase(i)
      dydx(6*nat+i) = pem(i)
      enddo
      do i=1,nsurft
      y(6*nat+nsurft+i) = cre(i)
      y(6*nat+2*nsurft+i) = cim(i)
      dydx(6*nat+nsurft+i) = gcre(i)
      dydx(6*nat+2*nsurft+i) = gcim(i)
      enddo

      if (methflag.eq.4) then
c     CSDM (coherent terms)
      do i=1,nsurft
      i2 = i + nsurft
      y(6*nat+3*nsurft+i) = cre(i2)
      y(6*nat+4*nsurft+i) = cim(i2)
      dydx(6*nat+3*nsurft+i) = gcre(i2)
      dydx(6*nat+4*nsurft+i) = gcim(i2)
      enddo
      endif

      endif

      return

      end
