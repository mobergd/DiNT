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

      
      subroutine initmpi
        
      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'mpif.h'

!      integer my_rank,nproc,ierr
      character(len=3) :: my_rank_str
      character(len=10) :: ldeltmp_str
      character(len=255) :: tmp_path

#if defined(MPITRAJS) || defined(MPIFORCES)
c     initialize mpi
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierr)
      write(my_rank_str,'(i0.3)') my_rank
#else
      my_rank = 0
      nproc = 1
#endif

#ifdef MPITRAJS
c     get the job_path
      call getcwd(job_path)

c     read tmp_path through DINT_TMP_DIR
      call getenv('DINT_TMP_DIR', tmp_path)
      if (tmp_path.eq.'') then
        ltmpdir = .false.
      else
        ltmpdir = .true.
      endif

c     read ldeltmp through DINT_DELETE_DIR
      call getenv('DINT_DELETE_TMP', ldeltmp_str)
      if (trim(ldeltmp_str).eq.'yes') then
        ldeltmp = .true.
      else
        ldeltmp = .false.
      endif

c     create work_path up to if DINT_TMP_DIR is set
      if (ltmpdir) then 
        call system( "mkdir "//trim(tmp_path)//" &> /dev/null" )
        work_path = trim(tmp_path)//"/"//my_rank_str 
        call system( "mkdir "//trim(work_path)//" &> /dev/null" )
      else
        work_path = trim(job_path)//"/"//my_rank_str 
        call system( "mkdir "//trim(work_path)//" &> /dev/null" )
      endif 

c     copy all files under job_path into work_path
      call system( "cp "//trim(job_path)//"/* "//trim(work_path)
     &//" &> /dev/null" )

c     change direction to work folder
      call chdir(work_path)
      call getcwd(tmp_path)
      if(trim(tmp_path).eq.trim(job_path)) then ! wrongly create folder, exit
        if (my_rank.eq.0) then
          write(*,*)
          write(*,*) "Wrong of enviroment variable $DINT_TMP_DIR."
          write(*,*) "E.g. if you set DINT_TMP_DIR=/home/work/tmp,"
          write(*,*) "make sure that upper folder /home/work exists."
          write(*,*)
        endif
        call MPI_BARRIER(MPI_COMM_WORLD, ierr)
        stop
      endif

c     open input file
!      open(5,FILE='input',ACTION='read')

c     open output file in path
!      open(6,FILE='output',ACTION='write')  
#endif

      end
