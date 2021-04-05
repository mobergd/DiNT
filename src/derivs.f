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

      
      subroutine derivs(x,yn,yout,nv,nsurf)

c This subroutine (along with XPTOY) puts together an single array 
c and its derivatives for integration.  It is called by the
c integrators whenever they need derivative information.  YN
c is the array of variables to be integrated, and YOUT is the
c derivatives.

      implicit none
      include 'param.f'
      include 'c_sys.f'

      integer nsurf,nv
      double precision xx(3,mnat),pp(3,mnat),gv(3,mnat),v,
     & cre(2*mnsurf),cim(2*mnsurf),gcre(2*mnsurf),gcim(2*mnsurf),
     & pem(mnsurf),phase(mnsurf),dvec(3,mnat,mnsurf,mnsurf)
      double precision yn(mnyarray),yout(mnyarray),phop(mnsurf)
      double precision x,dmag,gpem(3,mnat,mnsurf)

c     get xx from yn
      call xptoy(xx,pp,cre,cim,gv,gcre,gcim,mm,yn,yout,nat,nsurft,
     & 1,methflag,pem,phase)
c     get gradients
      call getgrad(xx,pp,nsurf,v,cre,cim,gv,gcre,gcim,nat,phop,dmag,
     & dvec,pem,gpem,phase)
c     transform to 1-D
      call xptoy(xx,pp,cre,cim,gv,gcre,gcim,mm,yn,yout,nat,nsurft,
     & 0,methflag,pem,phase)

      return
      end
