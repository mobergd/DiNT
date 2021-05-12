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


      subroutine diatom(jj,tau,rin,rout,mmm,xx,pp)

c     monte carlo selection of initial conditions for 
c     the atom-diatom initial conditions (INITx = 3)
c     Some of these routines were borrowed (with modifications)
c     from NAT8.1.

      implicit none
      include 'param.f'
      include 'c_ran.f'
#include <sprng_f.h>
c input
      double precision rin,rout,tau,jj
      double precision mmm(2),mi,mj,xx(3,2),pp(3,2)

c local
      integer i,irel,ii,ij
      real ran(5)
      double precision tin,tout,dist,rr,eta,cmom,thetap,phip,cosn
      double precision thetar,t,arg,xvec(5),tmp(12),r(3),r1(3),r2(3),
     & r3(3),rd(3)
     
      double precision phi,snt,cst,snp,csp,sne,cse,bmax,bmin,bsq
      double precision theta,mdum,rmass,xji,ppo,xxo

      do i=1,5
      xvec(i) = sprng(rng_stream)
      enddo

c one semiclassical choice
c      xji = jj+0.5d0
c a different choice
c      xji = dsqrt(jj*(jj+1.d0))
c yet a different choice
      xji = jj
c this code needs to be made consistent with respect to these choices

        xxo=2.0d0*xvec(1)-1.0d0
        phi=2.0d0*pi*xvec(2)
        eta=2.0d0*pi*xvec(3)
        t = tau*xvec(5)

        theta=dacos(xxo)
        tin = t
        tout = t-tau/2.0d0                                              2/4YL93
        if (tout.ge.0.0d0) go to 10                                     2/4YL93
        rr = rin
        go to 20
   10   continue
        rr = rout
   20 continue

c current version of the code starts the diatom at one turning point or the other...
c need to add sampling over phase

        snt=dsin(theta)
        cst=dcos(theta)
        snp=dsin(phi)
        csp=dcos(phi)
        sne=dsin(eta)
        cse=dcos(eta)
c ...   initialize coordinates
        tmp(1)=rr*snt*csp          
        tmp(2)=rr*snt*snp
        tmp(3)=rr*cst
c ...   initialize momenta
        ppo=xji/rr
        tmp(7)=-ppo*(snp*cse+csp*cst*sne)
        tmp(8)=ppo*(csp*cse-snp*cst*sne)
        tmp(9)=ppo*snt*sne

c transform Jacobis to x,y,z
        ii=1
        ij=2
        mi = mmm(ij)/(mmm(ii)+mmm(ij))
        mj = mmm(ii)/(mmm(ii)+mmm(ij))
        xx(1,ii) = -tmp(1)*mi
        xx(2,ii) = -tmp(2)*mi
        xx(3,ii) = -tmp(3)*mi
        xx(1,ij) = tmp(1)*mj
        xx(2,ij) = tmp(2)*mj
        xx(3,ij) = tmp(3)*mj

        pp(1,ii) = -tmp(7)
        pp(2,ii) = -tmp(8)
        pp(3,ii) = -tmp(9)
        pp(1,ij) = tmp(7)
        pp(2,ij) = tmp(8)
        pp(3,ij) = tmp(9)

        return

        end
