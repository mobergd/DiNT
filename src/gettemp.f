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


      subroutine gettemp(pp,mm,natom,temp,ke)

c Compute the temperature and kinetic energy from the momentum.

      implicit none
      include 'param.f'

      integer natom
      double precision pp(3,mnat),mm(mnat),temp,ke

c local
      integer i,j,k

      temp=0.d0
      ke=0.d0
      do i=1,natom
      do j=1,3
         ke = ke + 0.5d0*pp(j,i)*pp(j,i)/mm(i)
      enddo
      enddo
c <KE> = KE/Natom = 3/2 k T
      temp = 2.d0*ke/(kb*3.d0*dble(natom))
      return
      end
