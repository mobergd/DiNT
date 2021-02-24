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

      subroutine preptraj

c Put center of mass at origin

      implicit none
      include 'param.f'
      include 'c_traj.f'
      include 'c_sys.f'
#ifdef MPIFORCES
      include 'mpif.h'
#endif
!      integer my_rank,nproc,ierr
!      integer status(MPI_STATUS_SIZE)

      double precision xtot,ytot,ztot,mtot

      integer i

ccccc MPI
!      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
!      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
#ifdef MPIFORCES
      call MPI_BCAST(pp, 3*nat, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
      call gettemp(pp,mm,nat,temp,ke)
      IF (my_rank.eq.0) THEN
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF

      IF (my_rank.eq.0) THEN
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,nat
        xtot = xtot + mm(i)*xx(1,i)
        ytot = ytot + mm(i)*xx(2,i)
        ztot = ztot + mm(i)*xx(3,i)
        mtot = mtot + mm(i)
      enddo
      do i=1,nat
        xx(1,i) = xx(1,i) - xtot/mtot
        xx(2,i) = xx(2,i) - ytot/mtot
        xx(3,i) = xx(3,i) - ztot/mtot
      enddo
      ENDIF
#ifdef MPIFORCES
      call MPI_BCAST(xx, 3*nat, MPI_DOUBLE_PRECISION,
     &               0, MPI_COMM_WORLD, ierr)
      call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif

      call gettemp(pp,mm,nat,temp,ke)
      IF (my_rank.eq.0) THEN
      write(6,*)"Placed CoM at origin"
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF

c     remove overall momenta
      xtot=0.d0
      ytot=0.d0
      ztot=0.d0
      mtot=0.d0
      do i=1,nat
        xtot=xtot+pp(1,i)
        ytot=ytot+pp(2,i)
        ztot=ztot+pp(3,i)
        mtot = mtot + mm(i)
      enddo
      do i=1,nat
        pp(1,i)=pp(1,i)-xtot*mm(i)/mtot
        pp(2,i)=pp(2,i)-ytot*mm(i)/mtot
        pp(3,i)=pp(3,i)-ztot*mm(i)/mtot
      enddo
      call gettemp(pp,mm,nat,temp,ke)
      IF (my_rank.eq.0) THEN
      write(6,*)"Removed CoM motion"
      write(6,*)"Calculated temp = ",temp," K"
      write(6,*)"Total KE        = ",ke*autoev," eV"
      write(6,*)
      ENDIF

      return

      end
