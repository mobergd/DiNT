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


c COMMON BLOCKS FOR FREQUENTLY USED VARIABLES THAT ARE THE SAME FOR ALL
c TRAJECTORIES

c Note: When "i" is an index, the array runs over atom groups,
c   and when "j" is an index, the array runs over atoms,
c   and when "k" is an index, the array runs over surfaces.


c **********************************************************************
c  MPI VARIABLES
c **********************************************************************
c my_rank    Rank of processor
c my_itraj   The start # of the trajectory of this processor 
c my_ntraj   The end # of the trajectory of this processor
c nproc      Total number of prcoessor
c job_path   The absolute path of the folder you submit the job. This wll be got from pwd()
c work_path  The current working path for each processor. If DINT_TMP_DIR is defined, 
c            work_path equals to $DINT_TMP_DIR plus a processor number
c            e.g if DINT_TMP_DIR="/home/work", then for rank 1 processor,
c            there will be work_path="/home/work/001". If DINT_TMP_DIR is not defined,
c            work_path will be a created folder named by the rank number under the job_path.
c ltmpdir    If true, it means the DINT_TMP_DIR is not a null string 
c ldeltmp    If true, the tmp_path will be deleted after the calculation
c            This variable can be set by the sysytem enviroment variable 
c            DINT_DELETE_TMP. It can be set to be yes or no in lower letter.
c            e.g. " export DINT_DELETE_TMP=yes"  or
c                 " export DINT_DELETE_TMP=no"
        integer my_rank, my_itraj, my_ntraj, nproc
        character(len=256) :: work_path, job_path
        logical  ltmpdir,ldeltmp 
        common/c_mpi/my_rank, my_itraj, my_ntraj, nproc, work_path, 
     &          job_path, ltmpdir, ldeltmp



c **********************************************************************
c CONTROL VARIABLES
c **********************************************************************
c potflag    Flag that controls which potential interface to use
c nprint     Print info to unit 6 every nprint steps
c ntraj      Number of trajectories
c methflag   Flag for which single surface or non-BO method is used
c tflag      Special trajectory options
c repflag    Representation flag, =0 for adiabatic, =1 for diabatic
c ldofrag    If true, compute fragment energies. This option is useful 
c            for direct dynamics calculations where the fragments 
c            require different input options than the supermolecule.
c            Some options cannot be performed with this turned off.
c stodecoflag If true, use the stochastic decoherence.

      integer potflag,nprint,ntraj,methflag,tflag(4),repflag
      logical ldofrag,stodecoflag
      common/c_control/potflag,nprint,ntraj,methflag,tflag,repflag,
     &        ldofrag,stodecoflag
c **********************************************************************


c **********************************************************************
c SYSTEM
c **********************************************************************
c mm(j)       mass of atom j (input in amu, converted immediately to au)
c mmag(i)     total mass of atom group i (au, computed in readin)
c nsurft      number of electronic surfaces
c nmol        # of atom groups
c nat         # of total atoms (for all atom groups)
c natom(i)    # of atoms in AG i
c iatom(i)    index of first atom in AG i
c symbol(j)   atomic symbol of atom j
c ezero       zero of energy for the full potential (all atom groups), 
c             this number gets subtracted from the energies of every 
c             surface for all potential calls in GETPEM (input in eV, 
c             converted to au)
c ezeroim(i)  zero of energy for each AG, this value is used to correct
c             the zero of energy for each atom group (input in eV,
c             converted to au)

      double precision ezero
      double precision mm(mnat),mmag(mnmol),ezeroim(mnmol)
      integer nsurft,nmol,nat,natom(mnmol),iatom(mnmol)
      character*2 symbol(mnat)
      common/c_sys/ezero,ezeroim,mm,mmag,nsurft,nmol,nat,
     &             natom,iatom,symbol
c **********************************************************************


c **********************************************************************
c INITIAL CONDITIONS
c **********************************************************************
c nsurf0     Initial electronic surf for all trajectories
c initx(i)   Initial conditions flag for AG i (coords)
c initp(i)   Initial conditions flag for AG i (momenta)
c initj(i)   Initial conditions flag for AG i (ang mom)
c xx0(3,j)   Initial coordinates, read by readin (input in A, converted 
c            immediately to bohr)
c pp0(3,j)   Initial momenta, read by readin (input in au) 
c rran(3)    Numbers used to generate randum clusters if INITx = 1 (input in 
c            angstroms, converted immediately to bohr)
c temp0im(i) Target temperature of AG i when INITp = 0 (Kelvin)
c lescale0(i) If TRUE -> scale momenta to a target kinetic energy
c escale0im(i) Target kinetic energy of AG i when INITp = 0 and lescale0=TRUE 
c            (input in eV, immediately converted to hartree)
c scale0im(i) Target total energy of AG i when INITx = 5 and scale0im > 0

c The following are used for atom-diatom scattering conditions (INITx = 3)
c vvad        Initial vibrational quantum number of the diatom
c jjad        Initial rotational quantum number of the diatom
c arrad      Initial molecular arrangement
c escatad    Fixed total energy (input in eV, converted immediately to au)
c rrad        Initial atom-diatom distance (input in A, converted 
c            immediately to bohr)
c rinad      Inner turning point of initial diatom
c routad     Outer turning point of initial diatom
c tauad      Period in au of initial diatom
c bmaxad     Maximum impact parameter in au
c bminad     Minimum impact parameter in au
c pprelad    Initial relative momenta in au
c ecolad     Initial energy in au
c easym(k,m) Minimum energy for the diatom in arragement m
c rasym(k,m) Minimum energy distance for the diatom in arragement m

c The following are used for the methods based on normal modes
c nmqn(n,i)    Initial normal mode quanutm number for mode n
c nmtype(i)    Initial geometry is a well (0) or a saddle point (1).
c nmvec(3,j,m,i) Normal mode vectors for the m^th normal mode
c freq(m,i)    Normal mode frequencies for the m^th normal mode
c rturn(m,i)   Normal mode turning points for the m^th normal mode
c ewell(i)    Energy of the initial structure when using normal mode method for initial
c            conditions.  Should correspond to the energy of the bottom of some well.
c lreadhess  T=read Hessian from fort.70. F=compute Hessian numerically from gradients

c The following are used for diatom initial conditions (INITx = 4)
c vvdi(i)        Initial vibrational quantum number of the diatom for AG i
c jjdi(i)        Initial rotational quantum number of the diatom for AG i
c rmindi(i)      Guess at min diatomic separation for AG i, will be refined by code

c The following are used for reading initial conditions from a separate file (INITx=6)
c samptot(i)     Total number of samples contained in the separate file
c lbinsamp(i)    Samples files binary (direct access) or not
c sampfile(i)    Filename containing sampled data
c lems(i)        If true use EMS
c emicr(i)       Total energy of the microcanonical ensemble 
c nemstot(i)     Total number of points to output using EMS (efficient microcan. sampl.)
c ninc(i)        Incubation number of steps for EMS
c nbrea(i)       Number of steps in between points(trajectories) sampled in EMS
c emsw(i)        last calculated weight in EMS 
c lnorot(i)      If true, remove rotation. Equals to J=0 ensemble. 

c The following are used for photoexcited trajectories (TFLAG(3) = 1)
c ntarget        The target electronic state
c ephoton        The photon energy in au
c wphoton        The photon width / 2 in au (i.e., trajectories are excited 
c                if the electronic energy gap is ephoton +/- wphoton)

      integer nsurf0,initx(mnmol),initp(mnmol),initj(mnmol)
      double precision xx0(3,mnat),pp0(3,mnat),rran(3),
     &    temp0im(mnmol),escale0im(mnmol),scale0im(mnmol),emicr(mnmol)

      integer vvad,jjad,arrad
      double precision escatad,rrad,rinad,routad,
     & tauad,bmaxad,bminad,pprelad,ecolad,
     & easym(mnsurf,3),rasym(mnsurf,3)

      integer nmtype(mnmol)
      double precision nmvec(3,mnat,3*mnat,mnmol),freq(3*mnat,mnmol),
     & rturn(3*mnat,mnmol),ewell(mnmol),nmqn(3*mnat,mnmol)
      logical lreadhess,lbinsamp(mnmol),lems(mnmol),
     & letot(mnmol),lnorot(mnmol) 

      double precision rmindi(mnmol),vvdi(mnmol),jjdi(mnmol),
     &    rindi(mnmol),routdi(mnmol),taudi(mnmol),xtaudi(mnmol)

      integer samptot(mnmol),nemstot(mnmol),
     & ninc(mnmol),nbrea(mnmol)
      character*10 sampfilexx(mnmol),sampfilepp(mnmol)
      double precision samptarg(mnmol),sampjmin(mnmol),sampjmax(mnmol),
     & sampjtemp1(mnmol),sampjtemp2(mnmol),
     & sampjbrot1(mnmol),sampjbrot2(mnmol),
     & ejsc(4,mnmol),emsw(mnmol)

      integer ntarget
      double precision ephoton,wphoton

      common/c_initial/xx0,pp0,rran,
     & temp0im,escale0im,scale0im,emicr,escatad,rrad,rinad,routad,
     & tauad,bmaxad,bminad,pprelad,ecolad,easym,rasym,
     & nmvec,nmqn,freq,rturn,ewell,
     & rmindi,vvdi,jjdi,rindi,routdi,taudi,xtaudi,
     & ephoton,wphoton,
     & samptarg,sampjmin,sampjmax,sampjtemp1,sampjtemp2,emsw,
     & sampjbrot1,sampjbrot2,ejsc,
     & sampfilexx,sampfilepp,
     & ntarget,
     & nsurf0,initx,initp,initj,
     & vvad,jjad,arrad,
     & nmtype,
     & samptot,nemstot,ninc,nbrea,
     & lreadhess,lbinsamp,lems,letot
c **********************************************************************


c **********************************************************************
c INTEGRATOR VARIABLES
c **********************************************************************
c hstep0 = initial stepsize (input in fs, converted immediately to au)
c eps = converge energy for each step to this tolerance when using
c       the BS integrator (au)
c intflag = flag that controls which integrator to use

      double precision hstep0,eps
      integer intflag
      common/c_integrator/hstep0,eps,
     & intflag
c **********************************************************************

c **********************************************************************
c ORIENTATIONAL VARIABLES
c **********************************************************************
c iorient = determines how to orient AGs
c bmaxqc = maximium impact parameter
c bminqc = minimium impact parameter
c erelqc = relative energy
c tempqc = temperature of the relative energy distribution
c rel0qc = initial separation
c The following pertain to orienting the AGs with respect to each other
c comxx(i)   Initial CoM coords for AG i (input in A, converted
c            immediately to bohr)
c compp(i)   Initial CoM momenta for AG i (input in au)
      double precision
     & comxx(3,mnmol),compp(3,mnmol)
      integer iorient
      double precision bmaxqc,erelqc,rel0qc,bminqc,tempqc
      common/c_orient/comxx,compp,bminqc,bmaxqc,erelqc,tempqc,
     & rel0qc,iorient

c **********************************************************************
c TERMINATION CONDITIONS
c **********************************************************************
c t_stime = For TERMFLAG =1, Terminate after T_STIME fs.
c t_gradmag = For TERMFLAG =2, Terminate when the magnitude of the
c             gradient is less than T_GRADMAG in eV/A.
c t_r(o) = For TERMFLAG = 3, Termination distance, input in A, 
c          immediately converted to bohr.  
c t_nstep = Maximum number of steps for each trajectory.
c lwell = 0 geoms in one PE well are not selected, 1-acetylene, 2-vinyledene
c lchkdis = If true, check whether distances are less than 95% of the dissociation
c           termination distances. If the monitored distances are less than 95% of
c           the dissociation termination distances, the termination condition will 
c           not be checked and program will contiune to run.  

      double precision t_stime,t_gradmag,t_r(mnoutcome)
      integer termflag,t_nstep,t_noutcome,lwell
      character*2 t_symb(mnoutcome,2)
      logical lchkdis
      common/c_term/
     & t_r,t_stime,t_gradmag,termflag,t_nstep,t_noutcome,t_symb,lwell,
     & lchkdis
c **********************************************************************


c **********************************************************************
c SPECIAL TRAJECTORY OPTIONS
c **********************************************************************
c ramptime = ramp trajectory momentum every RAMPTIME time units 
c            (for TFLAG(1) = 2).  Read in fs, converted to au immediately.
c rampfact = ramp trajectory momentum by RAMPFACT (for TFLAG(1) = 2).
c nramp = ramp trajectory momentum NRAMP times (for TFLAG(1) = 2).
c andersen_temp = used by the Andersen thermostat, temperature
c andersen_freq = used by the Andersen thermostat, collision frequency
c scandth = Total energy to scale momenta after Andersen resampling (if >0)
      double precision ramptime,rampfact,andersen_temp,andersen_freq,
     & scandth
      integer nramp
      common/c_special/ramptime,rampfact,
     & andersen_temp,andersen_freq,scandth,
     & nramp
c **********************************************************************

c **********************************************************************
c MONTE CARLO INTEGRATION VARIABLES
c **********************************************************************
c mccurv = number of normal modes to treat as something other than a 
c          cartesian displacement (e.g., as a torsion or bend)
c mcmode = list of modes to be treated as a torsion, bend, etc.
c mctype = type of mode (1 = torsion)
c mcpar = additional parameters
      integer mccurv,mcmode(3*mnat),mctype(3*mnat),mcpar(2,3*mnat)
      common/c_montecarlo/mccurv,mcmode,mctype,mcpar
c **********************************************************************
