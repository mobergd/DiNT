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


c This file stores constants used for dimensioning arrays
c and for unit conversions.

c     max number of atom groups (aka molecules)
      integer mnmol
      parameter(mnmol=2)

c     max number of total atoms (total for all molecules)
      integer mnat
      parameter(mnat=100)

c     max number of electronic states
      integer mnsurf
      parameter(mnsurf=5)

c     max number of trajectories
      integer mntraj
      parameter(mntraj=1000000)

c     max number of termination condition outcomes
      integer mnoutcome
      parameter(mnoutcome=5)

c     maximum dimension of the Y-array that is integrated
      integer mnyarray
      parameter(mnyarray=6*mnat+5*mnsurf)

c     unit conversions
      double precision amutoau,kb,autoang,autoev,autofs,mu,autocmi,pi,
     & autokcal
      parameter(amutoau=1822.844987d0) ! best
c      parameter(amutoau=1822.888506d0) ! used by NAT
      parameter(kb=3.166829d-6)  ! Boltzmann in hartee/K
      parameter(autoang=0.52917706d0)
      parameter(autoev=27.2113961d0)
      parameter(autofs=0.024189d0)
      parameter(autocmi=219474.63067d0)
      parameter(autokcal=627.509d0)   ! IMLS
      parameter(mu=1.d0*amutoau)  ! mass-scaling mass
      parameter(pi=3.1415926536d0)
