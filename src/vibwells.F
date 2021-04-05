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


      subroutine vibwells(xxmt,ppmt,sceng,imm)

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_traj.f'
      include 'c_ran.f'
      include 'c_output.f'
#include <sprng_f.h>

      integer i,j,ii,imm,nfq
      double precision xxmt(3,mnat),ppmt(3,mnat),evb,sceng,
     & r1(mnat*3),tmp

      if (temp0im(imm).gt.0.d0) then
c       random thermal
        write(6,*)"Picking initial vibrational states",
     &   " for ",temp0im(imm)," K"
        write(6,*)"Mode #   Quanta"
        nfq=3*natom(imm)-6
        if (natom(imm).eq.2) nfq=1  ! diatom
c 201    evb=0.d0
        do i=1,natom(imm)
          ii = iatom(imm)+i
          do j=1,3
            xxmt(j,i) = xx0(j,ii)
            ppmt(j,i) = 0.d0
          enddo
        enddo
        do i=1,nfq
          call ranstate(freq(i,imm),nmqn(i,imm),temp0im(imm))
          write(6,*)i,nmqn(i,imm)
c          evb=evb+freq(i,imm)*(0.5d0+nmqn(i,imm))
        enddo
        write(6,*)
      else
c       random microcanonical
        nfq=3*natom(imm)-6
        do i=1,natom(imm)
          ii = iatom(imm)+i
          do j=1,3
            xxmt(j,i) = xx0(j,ii)
            ppmt(j,i) = 0.d0
          enddo
        enddo
        tmp = 0.d0
        do i=1,nfq
          r1(i) = sprng(rng_stream)
          tmp = tmp + r1(i)
        enddo
        write(6,*)"Mode QN and/or energy in unbound mode (au)"
        do i=1,nfq
          r1(i)=r1(i)/tmp*sceng   ! randomly distrubuted harmonic energy in mode i
          if (nmtype(imm).eq.1.and.i.eq.nfq) then
          nmqn(i,imm)=r1(i) ! Energy in unbound mode
          else
          nmqn(i,imm)=r1(i)/freq(i,imm)-0.5d0 ! effective quantum number for that much HO energy
          endif
          write(6,*)i,nmqn(i,imm)
        enddo
        write(6,*)"For a total harmonic energy of ",
     &  sceng*autoev," eV"
      endif
      return
          
      end
