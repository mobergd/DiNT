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

      subroutine jointprint

c This subroutine is to gather output of all processors' into a one
c The current path now should be located in the job path

      implicit none
      include 'param.f'
      include 'c_sys.f'
      include 'c_output.f'

      integer i,j,k,rank,io,fwrite,fread
      character(len=3)rank_str
      character(len=256) line,filep
      character(len=600) longline
      logical  existcode

ccccccccccccccccccccccccccccccccccccccccccc
c     unit 6 output: standard output
ccccccccccccccccccccccccccccccccccccccccccc
      if (.true.) then
        fwrite = 100 + 6
        fread  = 200 + 6 
        open(fwrite,FILE='output',ACTION='write')

c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/output'

c         open filep
          open(fread, FILE = filep, ACTION = 'read')        

c         print headstr of rank 0 file, ignore for other ranks
          read(fread,'(A)',iostat=io) line 
          
          do while(index(line," ***TRAJECTORY   ").eq.0)
c           check whether read flawlessly
            if (io.ne.0) then
              write(fwrite, *) "Wrong when reading file ",filep
              return
            endif
c           only for the rank 0
            if (rank.eq.0) then 
c              print ranks' corresponding trajectories number
              if (index(line,"...entering READIN...").ne.0) then
                write( fwrite, * )''
                write(fwrite,'("... MPI environment setting ...",/)')
                do i=0,nproc-1
                  write(fwrite,'(" rank ",i3," ran trajectories",
     & i6," to",i6)')i,i*int(ceiling(ntraj/dble(nproc)))+1,
     & min(ntraj, (i+1)*int(ceiling(ntraj/dble(nproc))))
                enddo
                write(fwrite,'(/," Job-submitted path = ",A,/)')
     &trim(job_path)
                write(fwrite,'(" Work path = ",A,/)') 
     &trim(work_path)

                if (ldeltmp) then
                  write(fwrite,'(" Delete the work  
     &folder after calcultion: Yes "/)')
                else
                  write(fwrite,'(" Delete the work 
     &folder after calculation: No "/)')
                endif
              endif
              write(fwrite,'(A)') trim(line)
            endif
            read(fread, '(A)', iostat=io) line
          enddo

c         the same print for all ranks
          do while(io.eq.0)
            write(fwrite,'(A)') trim(line)
            read(fread,'(A)',iostat=io) line
          enddo

        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 10 output
ccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(10) ) then
        fwrite = 100 + 10
        fread  = 200 + 10
        open(fwrite,FILE='fort.10',ACTION='write')
        
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.10'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then          
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           hederstr
            do i=0,1
              read( fread, '(A)', iostat = io) line
              if (rank.eq.0) then
                if (io.ne.0) then
                  write( fwrite, * ) "Wrong when reading file ",filep
                  return
                endif
                write(fwrite,'(A)') trim(line)
              endif
            enddo

c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 11 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(11) ) then
        fwrite = 100 + 11
        fread  = 200 + 11
        open(fwrite,FILE='fort.11',ACTION='write')

c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.11'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then  
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )
c           hederstr
            do i=0,1
              read( fread, '(A)', iostat = io) line
              if (rank.eq.0) then
                if (io.ne.0) then
                  write( fwrite, * ) "Wrong when reading file ",filep
                  return
                endif
                write(fwrite,'(A)') trim(line)
              endif
            enddo      
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 20 output
ccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(20) ) then
        fwrite = 100 + 20
        fread  = 200 + 20
        open(fwrite,FILE='fort.20',ACTION='write')

c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.20'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then            
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           hederstr
            do i=0,1
              read( fread, '(A)', iostat = io) line
              if (rank.eq.0) then
                if (io.ne.0) then
                  write( fwrite, * ) "Wrong when reading file ",filep
                  return
                endif
                write(fwrite,'(A)') trim(line)
              endif
            enddo
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 21 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(21) ) then
        fwrite = 100 + 21
        fread  = 200 + 21
        open(fwrite,FILE='fort.21',ACTION='write')
      
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.21'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then          
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           hederstr
            do i=0,1
              read( fread, '(A)', iostat = io) line
              if (rank.eq.0) then
                if (io.ne.0) then
                  write( fwrite, * ) "Wrong when reading file ",filep
                  return
                endif
                write(fwrite,'(A)') trim(line)
              endif
            enddo
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 30 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(30) ) then
        fwrite = 100 + 30
        fread  = 200 + 30
        open(fwrite,FILE='fort.30',ACTION='write')
      
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.30'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then          
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           the same print for all ranks
            read(fread,'(A)',iostat=io) longline
            do while(io.eq.0)
              write(fwrite,'(A)') trim(longline)
              read(fread,'(A)',iostat=io) longline
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 31 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(31) ) then
        fwrite = 100 + 31
        fread  = 200 + 31
        open(fwrite,FILE='fort.31',ACTION='write')
      
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.31'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then          
c           open filep
            open(fread, FILE = filep, ACTION = 'read')
c           the same print for all ranks
            read(fread,'(A)',iostat=io) longline
            do while(io.eq.0)
              write(fwrite,'(A)') trim(longline)
              read(fread,'(A)',iostat=io) longline
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 40 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(40) ) then
        fwrite = 100 + 40
        fread  = 200 + 40
        open(fwrite,FILE='fort.40',ACTION='write')
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.40'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then
c           open filep
            open(fread, FILE = filep, ACTION = 'read')
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 41 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(41) ) then
        fwrite = 100 + 41
        fread  = 200 + 41
        open(fwrite,FILE='fort.41',ACTION='write')
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.41'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then  
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           the same print for all ranks
            do while(io.eq.0)
              read(fread,'(A)',iostat=io) line
              write(fwrite,*) trim(line)
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 42 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(42) ) then
        fwrite = 100 + 42
        fread  = 200 + 42
        open(fwrite,FILE='fort.42',ACTION='write')
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.42'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then  
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 43 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(43) ) then
        fwrite = 100 + 43
        fread  = 200 + 43
        open(fwrite,FILE='fort.43',ACTION='write')
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.43'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then  
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
cccccccccccccccccccccccccccccccccccccccccc
c     unit 80 output
cccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(80) ) then
        fwrite = 100 + 80
        fread  = 200 + 80
        open(fwrite,FILE='fort.80',ACTION='write')
      
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.80'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then  
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           hederstr
            do i=0,1
              read( fread, '(A)', iostat = io) line
              if (rank.eq.0) then
                if (io.ne.0) then
                  write( fwrite, * ) "Wrong when reading file ",filep
                  return
                endif
                write(fwrite,'(A)') trim(line)
              endif
            enddo
          
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif
ccccccccccccccccccccccccccccccccccccccccccc
c     unit 81 output
ccccccccccccccccccccccccccccccccccccccccccc
      if ( lwrite(81) ) then
        fwrite = 100 + 81
        fread  = 200 + 81
        open(fwrite,FILE='fort.81',ACTION='write')
      
c       loop over all folders
        do rank=0,nproc-1
c         set processor-specific paramters     
          write(rank_str,'(i0.3)') rank
          filep = trim(job_path)//'/'//rank_str//'/fort.81'
          inquire(file=filep, exist=existcode)
          if (existcode.eqv..true.) then  
c           open filep
            open( fread, FILE = filep, ACTION = 'read' )      
c           hederstr
            do i=0,1
              read( fread, '(A)', iostat = io) line
              if (rank.eq.0) then
                if (io.ne.0) then
                  write( fwrite, * ) "Wrong when reading file ",filep
                  return
                endif
                write(fwrite,'(A)') trim(line)
              endif
            enddo
c           the same print for all ranks
            read(fread,'(A)',iostat=io) line
            do while(io.eq.0)
              write(fwrite,*) trim(line)
              read(fread,'(A)',iostat=io) line
            enddo
          endif
        enddo
      endif 

      end
