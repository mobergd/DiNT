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

      subroutine premol(im)

c Precompute information for each atom group.

      implicit none
      include 'param.f'
      include 'c_sys.f'
 
      double precision xxm(3,mnat),ppm(3,mnat),mmm(mnat)
      double precision nmvecm(3,mnat,3*mnat),freqm(3*mnat),
     & rturnm(3*mnat),ewellm,nmqnm(3*mnat)
      integer im,i,j,ii,k
      character*2 symb(mnat)

      write(6,*)"Precomputing info for AG ",im
      write(6,*)"------------------------------"

c     store coords for this molecule in their own array for easier manipulation
      do i=1,natom(im)
       ii=i+iatom(im)
       mmm(i) = mm(ii)
       symb(i) = symbol(ii)
       do j=1,3
        xxm(j,i) = 0.d0
        ppm(j,i) = 0.d0
       enddo
      enddo

      if (initx(im).eq.2.or.initx(im).eq.5) then
c       SET UP FOR NORMAL MODE METHODS
c       set to initial structure
        do i=1,natom(im)
        ii=i+iatom(im)
        mmm(i) = mm(ii)
        do j=1,3
          xxm(j,i) = xx0(j,ii)
        enddo
        enddo

        do i=1,natom(im)*3
          nmqnm(i) = nmqn(i,im)
        enddo
        call normod(xxm,mmm,natom(im),nmvecm,freqm,ewellm,
     &   repflag,nsurf0,rturnm,nmqnm,lreadhess,ldofrag,
     &   symb,nmtype(im))
        ewell(im) = ewellm
        do i=1,3*natom(im)
          freq(i,im)=freqm(i)
          rturn(i,im)=rturnm(i)
          do j=1,3
          do k=1,natom(im)
            nmvec(j,k,i,im)=nmvecm(j,k,i)
          enddo
          enddo
        enddo
      else if (initx(im).eq.3) then
c ATOM-DIATOM METHOD
        write(6,*)"INITx = 3:  Atom-diatom initial conditions"
        call preatomdiatom(escatad,jjad,vvad,arrad,
     &    rrad,mmm,rinad,routad,tauad,bmaxad,
     &    bminad,pprelad,ecolad,nsurf0,nsurft,easym,rasym)
      else if (initx(im).eq.4) then
c DIATOM METHOD
        write(6,*)"INITx = 4:  Diatom initial conditions"
        call prediatom(vvdi(im),jjdi(im),rmindi(im),mmm,nsurf0,
     &     symb,taudi(im),rindi(im),routdi(im))
      else
        write(6,*)"Nothing to do..."
        write(6,*)
      endif

      return

      end
