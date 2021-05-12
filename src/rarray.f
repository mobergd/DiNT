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


      subroutine rarray(ithistraj,time,xx,nat)

      implicit none
      include 'param.f'

      integer i,j,k,nat,l,ithistraj
      double precision rr(mnat*(mnat-1)/2),rij,xx(3,mnat),time,
     & symbol(mnat)

      l=0
      do 10 i=1,nat
      do 10 j=i+1,nat
      rij=0.d0
      do k=1,3
      rij = rij + (xx(k,i)-xx(k,j))**2
      enddo
      rij = dsqrt(rij)
      l=l+1
      rr(l)=rij
 10   continue

      write(43,143)ithistraj,time*autofs,(rr(k)*autoang,k=1,l)
 143  format(i10,f20.8,100(f20.8))

      return
      end
