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


      subroutine getrho(cre,cim,rhor,rhoi,nsurft)

c Compute the density matrix from the electronic coefficients.
c Note:  Phases are handled separately, and RHOR and RHOI do
c not include phases.  These may need to be added to RHOR and RHOI,
c depending on how these quantities are used.

      implicit none
      include 'param.f'

      integer nsurft,i,j
      double precision cim(2*mnsurf),cre(2*mnsurf),
     & rhor(mnsurf,mnsurf),rhoi(mnsurf,mnsurf)

      do i=1,nsurft
      do j=1,nsurft
        rhor(i,j) = cre(i)*cre(j)+cim(i)*cim(j)
        rhoi(i,j) = cre(i)*cim(j)-cre(j)*cim(i)
      enddo
      enddo

      return

      end
