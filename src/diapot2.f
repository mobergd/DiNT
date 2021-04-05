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


      subroutine diapot2(x,v,xj,mmm,nsurf,symb)

c Compute the diatomic energy.

      implicit none
      include 'param.f'
      include 'c_sys.f'
      integer j,nsurf,i
      double precision v,xx(3,mnat),xj,erot,rmass,x,
     &  mmm(2)

      double precision pemd(mnsurf,mnsurf),gpemd(3,mnat,mnsurf,mnsurf),
     & pema(mnsurf),gpema(3,mnat,mnsurf),dvec(3,mnat,mnsurf,mnsurf)

      character*2 symb(mnat)

      xx(1,1) = 0.d0
      xx(2,1) = 0.d0
      xx(3,1) = 0.d0
      xx(1,2) = 0.d0
      xx(2,2) = 0.d0
      xx(3,2) = x
      call getpem(xx,2,pema,pemd,gpema,gpemd,dvec,symb)
      if (repflag.eq.1) then
        v = pemd(nsurf,nsurf)
      else
        v = pema(nsurf)
      endif

c      write(6,*)x*autoang,v

c add rotation if j > 0
      if (xj.ge.0) then
       rmass = mmm(1)*mmm(2)/(mmm(1)+mmm(2))
       erot = 0.5d0*(xj+0.5d0)**2/(rmass*x**2)
       v = v + erot
      endif

      return
      end

