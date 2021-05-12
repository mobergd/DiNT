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

      subroutine endmpi
        implicit none
        include 'param.f'
        include 'c_sys.f'
        include 'mpif.h'

        integer err_code

        if (ltmpdir) then
c     copy folder work_path back to job_path
          call system( "cp -R "//trim(work_path)//"  "//trim(job_path) )
c     clear the tmp_path
          if (ldeltmp) then
            call system( "rm -rf "//trim(work_path) )
          endif
        endif

c     barrier
        call MPI_Barrier(MPI_COMM_WORLD, err_code)

c     joint output
        if (my_rank.eq.0) then
c     change direction to job folder
          call chdir(job_path)
          call jointprint
        endif

        call MPI_Finalize( err_code )
      end
