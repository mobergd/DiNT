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
      
      subroutine bootstrap(dd,wd,nd,md,id,perr)

      implicit none
c      include 'param.f'
c      include 'c_ran.f'
c#include <sprng_f.h>

      integer nd,md,id,iseed,nsamp,irand,j,i
      double precision dd(md),avgi(100000),x,avgavg,tot,wtot,
     & avg,std,perr,wd(md),rand
      integer*4 timeArray(3)

      call itime(timeArray)     ! Get the current time
      iseed = timeArray(1)+timeArray(2)+timeArray(3)
      x=rand(iseed)

      nsamp=1000

      avgavg=0.d0
      do j=1,nsamp
        tot=0.d0
        do i=1,id
          irand=int(rand()*dble(id))+1
          if (irand.le.nd) then
              tot=tot+dd(irand)*wd(irand)
          endif
        enddo
        avg = tot/dble(id)
        avgi(j)=avg
        avgavg=avgavg+avg
      enddo
      avgavg=avgavg/dble(nsamp)

      std=0.
      do j=1,nsamp
        std=std+(avgi(j)-avgavg)**2
      enddo
      std=std/dble(nsamp-1)
      std=dsqrt(std)

      perr=std/avg*100.d0
c      write(99,*)avg,std,std/avg*100.d0

      return
      end
