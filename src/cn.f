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


      subroutine cn(xx,nclu,cntot,time)

      implicit none
      include 'param.f'

      integer i,j,k,nclu
      double precision xx(3,mnat),lind,rij,cnx,del,rnn,rcut,
     &  cntot,time,rcom

      rnn = (4.0217/dsqrt(2.d0))/autoang
      rcut = rnn*1.5d0
      del = rnn/2.d0

      cntot = 0.d0
      do i=1,nclu
        cnx = 0.d0
c compute distance of atom i from center of mass (assumed to be at the origin)
c        rcom = 0.d0
c        do k=1,3
c        rcom = rcom + xx(k,i)**2
c        enddo
c        rcom = dsqrt(rcom)
c compute coordination number for atom i
        do j=1,nclu
          if (i.ne.j) then
          rij = 0.d0
          do k=1,3
          rij = rij + (xx(k,i)-xx(k,j))**2
          enddo
          rij = dsqrt(rij)
          cnx = cnx + (1.d0-dtanh(del*(rij-rcut)))/2.d0
          endif
        enddo
c     compute average CN
      cntot = cntot + cnx
      enddo
      cntot = cntot/nclu

      return
 
      end
