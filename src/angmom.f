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
      
      subroutine angmom(xx,pp,mm,natom,bigj,bigjtot)

c computes angular momentum around the origin using the position XX and momentum PP
c components of angmom are the cartesian components, not the components around the the principal axes

      implicit none
      include 'param.f'

      integer natom
      double precision xx(3,natom),pp(3,natom),mm(natom),
     &   bigj(3),bigjtot

c local
      integer i

c zero
      bigjtot=0.d0
      do i=1,3
       bigj(i)=0.d0
      enddo

c bigj = xx cross pp (J = R x P)
      do i=1,natom
c       write(6,*)"xx",(xx(k,i),k=1,3),"pp",(pp(l,i),l=1,3)
       bigj(1)=bigj(1)+(xx(2,i)*pp(3,i)-xx(3,i)*pp(2,i))
       bigj(2)=bigj(2)-(xx(1,i)*pp(3,i)-xx(3,i)*pp(1,i))
       bigj(3)=bigj(3)+(xx(1,i)*pp(2,i)-xx(2,i)*pp(1,i))
      enddo

c bigj(i) is the ith component of the total angular momentum
c bigjtot is the magnitide of the total angular momentum
      do i=1,3
       bigjtot = bigjtot + bigj(i)**2
      enddo
      bigjtot = dsqrt(bigjtot)

      end
