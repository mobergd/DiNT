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


      subroutine ewkb2(ji,vi,mm,erotvib,rin,rout,nsurf,symb,rguess)

c This subroutine is specialized for the atom-diatom initial conditions.
c It computes the WKB energy for the rovibrational state (vi,ji).
c INPUT
c ji = initial rot state of diatom
c vi = initial vib state of the diatom
c mm = masses of the atoms
c OUTPUT
c erotvib = rivibrational energy for (vi,ji) state

      implicit none
      include 'param.f'

      integer arri,nsurf
      double precision mm(mnat),erotvib,xj,ji,vi
      character*2 symb(mnat)

c     local
      integer i,j,JMAX
      double precision rmin,emin,emax,xv,f,fmid,dx,rtbis_vib,etry,
     & rin,rout,rguess,etmp

      xj = ji

c      write(6,*)"aj in ewkb2"

c Need to bracket energy for (vi,ji) state
c     find minimum of diatom with rotation
      call diamin2(rmin,emin,xj,mm,nsurf,symb,rguess)
c     get asymptotic energy
      call diapot2(30.d0,emax,0.d0,mm,nsurf,symb)
c      write(6,*)"Minimum energy for diatom ",
c     &  " is ",emin*autoev," eV at ",rmin*autoang," A",xj
c      write(6,*)"Asym energy is ",emax*autoev," eV"
c      write(6,*)

c Now determine energy of (vi,ji) state
c     compute vibrational action at E = emax
      rin=0.5d0
      rout=30.d0
      call vwkb2(mm,emin,xj,xv,rin,rout,nsurf,symb)
      f = xv - vi
c      write(6,*)"vibrational action at emin = ",xv
c     compute vibrational action at E = emin
      rin=0.5d0
      rout=30.d0
      etmp=0.99d0*emax
      call vwkb2(mm,etmp,xj,xv,rin,rout,nsurf,symb)
      fmid = xv - dble(vi)
c      write(6,*)"vibrational action at 0.99 emax = ",xv

      if (f*fmid.ge.0.d0) then
          write(6,1000)
          stop
      end if
      if (f.lt.0.d0) then
          rtbis_vib=emin
          dx=etmp-emin
      else
          rtbis_vib=etmp
          dx=emin-etmp
      end if

      JMAX=100
      do 11 j = 1, JMAX
          dx=dx*0.5d0
          etry=rtbis_vib+dx
          rin=0.5d0
          rout=30.d0
          call vwkb2(mm,etry,xj,xv,rin,rout,nsurf,symb)
c          write(6,*)"trying...",etry*27.211d0
c          write(6,*)"vibrational action at etry = ",xv,vi
          fmid = xv - vi
          if (fmid.le.0.d0) rtbis_vib=etry
          if (abs(dx).lt.1.d-10.or.fmid.eq.0.d0) go to 12
   11 continue
      write(6,1010)
      stop
 1000 format(2x,'root must be bracketed in rtbis_vib')
 1010 format(2x,'Too many bisections in rtbis_vib')

   12 continue

      erotvib = rtbis_vib

      return

      end
