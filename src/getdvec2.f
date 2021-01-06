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

      subroutine getdvec2(nclu,u2b2,gu2b2,dvec2b2)

c Computes an effective nonadiabatic coupling vector c
c from a 2x2 diabatic matrix.

      implicit none
      include 'param.f'

      integer i1,i2,k,l,nclu
      double precision u2b2(2,2),gu2b2(3,mnat,2,2),dvec2b2(3,mnat)
      double precision u11,u12,u22,u21,cc(2,2),v1,v2,sn1,cs1,tmp

c     compute adiabatic info
c     diagonalize the 2x2
      u11 = u2b2(1,1)
      u12 = u2b2(1,2)
      u21 = u2b2(2,1)
      u22 = u2b2(2,2)
      call dlaev2(u11,u12,u22,v2,v1,cs1,sn1)
c     phase convention (this part depends on how DLAEV2 works)
      if (v2.ge.v1) then
        cc(1,1) = -sn1
        cc(2,1) = cs1
        cc(1,2) = cs1
        cc(2,2) = sn1
      else
        tmp = v1
        v1 = v2
        v2 = tmp
        cc(1,1) = cs1
        cc(2,1) = sn1
        cc(1,2) = -sn1
        cc(2,2) = cs1
      endif

c compute d
      do i1 = 1,3
      do i2 = 1,nclu
        dvec2b2(i1,i2) = 0.d0
        do k = 1, 2
        do l = 1, 2
          dvec2b2(i1,i2) = dvec2b2(i1,i2)
     &                    +cc(k,2)*cc(l,1)*gu2b2(i1,i2,k,l)
        enddo
        enddo
        if ((v1-v2) .ne. 0.0d0) then
          dvec2b2(i1,i2) = dvec2b2(i1,i2)/(v2-v1)
        else
          dvec2b2(i1,i2) = 0.0d0
        endif
      enddo
      enddo

      return

      end
