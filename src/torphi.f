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


      subroutine torphi(xx,nclu,a1,a2,phi,n1)

      implicit none
      include 'param.f'

      integer nclu,a1,a2,n1
      double precision xx(3,mnat),phi

c local
      integer i,j
      double precision axis(3),xx2(3,mnat),cp,op,sp,tmp,rot(3,3),
     & tmp1,tmp2,tmp3

c move to put one of the atoms at the origin
      do i=1,nclu
      do j=1,3
        xx2(j,i)=xx(j,i)-xx(j,a2)
      enddo
      enddo

c rotation axis
      tmp=0.d0
      do i=1,3
        axis(i)=xx2(i,a1)-xx2(i,a2)
        tmp=tmp+axis(i)**2
      enddo
      do i=1,3
      axis(i)=axis(i)/dsqrt(tmp)
      enddo

c rotation matrix
c      cp=dcos(phi/180.d0*pi)
c      sp=dsin(phi/180.d0*pi)
      cp=dcos(phi)
      sp=dsin(phi)
      op=1.d0-cp

      rot(1,1)=cp+axis(1)**2*op
      rot(2,2)=cp+axis(2)**2*op
      rot(3,3)=cp+axis(3)**2*op
      rot(1,2)=axis(1)*axis(2)*op-axis(3)*sp
      rot(2,1)=axis(1)*axis(2)*op+axis(3)*sp
      rot(1,3)=axis(1)*axis(3)*op+axis(2)*sp
      rot(3,1)=axis(1)*axis(3)*op-axis(2)*sp
      rot(2,3)=axis(2)*axis(3)*op-axis(1)*sp
      rot(3,2)=axis(2)*axis(3)*op+axis(1)*sp

c rotate
      do i=1,n1
         tmp1 = xx2(1,i)
         tmp2 = xx2(2,i)
         tmp3 = xx2(3,i)
         xx2(1,i)=tmp1*rot(1,1)+tmp2*rot(2,1)+tmp3*rot(3,1)
         xx2(2,i)=tmp1*rot(1,2)+tmp2*rot(2,2)+tmp3*rot(3,2)
         xx2(3,i)=tmp1*rot(1,3)+tmp2*rot(2,3)+tmp3*rot(3,3)
      enddo

c move back 
      do i=1,nclu
      do j=1,3
        xx2(j,i)=xx2(j,i)+xx(j,a2)
      enddo
      enddo
      do j=1,3
      do i=1,nclu
        xx(j,i)=xx2(j,i)
      enddo
      enddo

      return

      end
