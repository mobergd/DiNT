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


      subroutine vwkb(arr,mm,edia,xj,xv,rin,rout,nsurf)

c Computes the WKB vibrational action at the energy EDIA and
c for the rotational state XJ.  For the atom-diatom initial
c conditions method (INITx = 3) only.  From NAT8.1.

      implicit none
      include 'param.f'

      integer arr,ncheb,i,nsurf
      double precision mm(mnat),edia,xv,sum,arg,x(12),wgt(12),
     &  v,rin,rout,rmass,term,r,rmid,rdif,xj
      
      ncheb = 12
      do 5 i = 1,ncheb
        arg = dble(i)/dble(ncheb+1)*pi
        x(i) = dcos(arg)
        wgt(i) = pi/dble(ncheb+1)*(dsin(arg))**2
    5 continue

      call turn(arr,mm,edia,xj,rin,rout,nsurf)
      if (rin .gt. rout) then                                  
        write(6,*)'Outer turning point smaller than inner turning',
     &      ' point in vwkb, trajectory'
        stop
      else if ((rout - rin) .lt. 1.d-5) then
        xv = -0.5d0                    
      else                              
        rdif = 0.5d0*(rout-rin)
        rmid = 0.5d0*(rout+rin)
        sum = 0.d0
        do 10 i = 1,ncheb
          r = rmid+rdif*x(i)
          call diapot(r,arr,v,xj,mm,nsurf)
          if( (edia-v)/((r-rin)*(rout-r)) .lt. 0.0d0)then
            write(6,*)'WARNING vwkb,  attempting ',
     &       'to take a sqrt of a negative number ',edia-v
            write(6,*)'Change guess in diamin?'
            stop
          else
            term = sqrt((edia-v)/((r-rin)*(rout-r)))
          endif
          sum = sum + wgt(i)*term
   10   continue
      if (arr.eq.1) rmass = mm(1)*mm(2)/(mm(1)+mm(2))
      if (arr.eq.2) rmass = mm(2)*mm(3)/(mm(2)+mm(3))
      if (arr.eq.3) rmass = mm(3)*mm(1)/(mm(3)+mm(1))
        xv = sqrt(2.d0*rmass)*rdif**2*sum/pi-0.5d0
      endif

      return

      end
