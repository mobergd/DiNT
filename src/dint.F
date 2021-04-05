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


      program dint
c     Dint:  Adiabatix and Nonadiabatic Trajectories
c     Main program.  This routine calls HEADER and READIN once, then
c     calls PREMOL once per AG to precalculate some constants.
c     DRIVER is called NTRAJ times, once for each trajectory.
c     FINALSTATE is called after each trajectory to write some data.

c     MPI version of Dint modified by Rui Ming Zhang


      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_ran.f'
      include 'c_traj.f'
      include 'mpif.h'
#include <sprng_f.h>

      integer im,initrng,minitrng,i,ia, n_tmp
      double precision t1,t2,t3,cputime,t3old,r

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     INITIALIZE DINT2.0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     init MPI
      call initmpi

c     print header
      call header

      write(6,'("...entering READIN...",/)')
      call readin
      write(6,'("...exiting READIN...",/)')

c     prepare the pot.
c     do nothing for analytical PES
c     read qc.txt for dd
      write(6, '("... entering PREOPT ... ",/)')
      call prepot
      write(6, '("... exiting PREOPT ... ",/)')

c     prepare each molecule
c     compute quantities that are the same for all trajectories
c     initialize electronic coordinates
      call timing(t1)
      write(6,'("...entering PREMOL...",/)')
      call initelec
      do im = 1 , nmol
      call premol(im)
      enddo
      write(6,'("...exiting PREMOL...",/)')
      call timing(t2)
      write(6,'(" CPU time in PREMOL is ",f10.5," s",/)')(t2-t1)

c     assign trajectories for each processor
      if (ntraj.lt.nproc) then
        write(*,*)" Wrong, the number of processor exceeds the 
     &number of trajectories."
        stop
      endif
      n_tmp    = int( ceiling( ntraj/dble(nproc) ) )
      my_itraj = my_rank * n_tmp + 1
      my_ntraj = min( ntraj, ( my_rank + 1 ) * n_tmp )
      
      t3 = 0.d0
      cputime = 0.d0

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     LOOP OVER TRAJECTORIES
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      nistep=0
      nistepw=0

      do itraj = my_itraj,my_ntraj

      write(6,*)"***TRAJECTORY ",itraj," OUT OF ",ntraj,"***"
      write(6,*)

      if (tflag(2).eq.1) then
      write(6,*)"Restarting trajectory #",trajlist(itraj)
      write(6,*)
      initrng = trajlist(itraj)-1
      else
      initrng = itraj-1
      endif

c     initialize this trajectory
      nsurf = nsurf0
      hstep = hstep0
c     initialize electronic coordinates
      call initelec
c     initialize the randum number generator
      rng_stream = init_sprng(initrng,maxtraj,ranseed,SPRNG_DEFAULT)

      write(6,'("...entering INITMOL...",/)')
      do im = 1 , nmol
c     generate specific initial coordinates for each molecule
      call initmol(im)
      enddo
      write(6,'("...exiting INITMOL...",/)')

      write(6,'("...entering ORIENT...",/)')
      do im = 1 , nmol
c     orient AGs with respect to each other
      call orient(im)
      enddo
      write(6,'("...exiting ORIENT...",/)')

      write(6,'("...entering PREPTRAJ...",/)')
      call preptraj
      write(6,'("...exiting PREPTRAJ...",/)')

c     propagate this trajectory
      write(6,'("...entering DRIVER...",/)')
      call driver
      write(6,'("...exiting DRIVER...",/)')
      t3old = t3
      call timing(t3)
      cputime = cputime + (t3-t3old)
      write(6,'(" CPU time for this trajectory is ",f10.5," s",/)')
     & (t3-t3old)

c     analyze trajectory
      write(6,'("...entering FINALSTATE...",/)')
      call finalstate
      write(6,'("...exiting FINALSTATE...",/)')

c     end loop over trajectories
      enddo

      write(6,'(" CPU time for all trajectories is ",f10.5," s or ",
     & f10.5," s per trajectory.",/)')
     & cputime,cputime/dble(ntraj)

 999  write(6,'("Well done!")')

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     END 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      call endmpi

      end
