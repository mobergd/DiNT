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


      subroutine prediatom(vi,ji,rguess,mm,nsurf,symb,
     &    tau,rin,rout)

c Precompute some information for the atom-diatom initial conditions
c (INITx = 4).

      implicit none
      include 'param.f'

      integer nsurf,arri
      double precision mm(mnat),escat,rguess,
     & einti,eroti,evibi,eelec,rmin,rin,rout,ecoli,rmass,term,
     & mdum,ppreli,bmax,bmin,tau,etmp,
     & ji,vi
      character*2 symb(mnat)

      write(6,*)"Asymptotic analysis"
      write(6,*)"   Emin (eV)      Rmin (A)"
      call diamin2(rmin,etmp,-1.d0,mm,nsurf,symb,rguess)
      write(6,50)etmp*autoev,rmin*autoang
 50   format(2f15.5)

      write(6,*)"Preparing the diatom ",
     &  "with v = ",vi," j = ",ji," surface = ",nsurf
      write(6,*)
      eelec = etmp
      write(6,*)"Energetics (relative to zero of E)"
      write(6,*)"At the classical minimum:"
      write(6,100)"Re    = ",rmin*autoang," A"
      write(6,100)"Eelec = ",eelec*autoev," eV"
      call ewkb2(0.d0,vi,mm,evibi,rin,rout,nsurf,symb,rmin)
      call ewkb2(ji,vi,mm,einti,rin,rout,nsurf,symb,rmin)
      eroti = einti-evibi

      write(6,102)"Eint  = E(v=",vi,",j=",ji,") =",einti*autoev," eV"
      write(6,102)"Evib  = E(v=",vi,",j=",0.,") =",evibi*autoev," eV"
c     & ,2.d0*evibi*autocmi," cm-1"
      write(6,103)"Erot  = Eint  - Evib   = ",eroti*autoev," eV"
      write(6,101)"Turning points = ",rin*autoang," and ",
     & rout*autoang," A"
 100  format(1x,a,20x,f13.5,a)
 101  format(1x,a,f13.5,a,f13.5,a)
 102  format(1x,a,f6.2,a,f6.2,a,4x,f13.5,a,f13.5,a)
 103  format(1x,a,3x,f13.5,a)

      call period2(einti,ji,mm,rin,rout,tau,nsurf,symb)
      write(6,103)"Period for the vibration = ",tau*autofs," fs"
      write(6,*)

      return

      end
