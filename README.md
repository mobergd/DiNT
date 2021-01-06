# DiNT
## Direct Nonadiabatic Trajectories v2.0: A code for non-Born--Oppenheimer molecular dynamics

This is the development version of DiNT obtained from the AutoMech suite: [AutoMech](https://github.com/Auto-Mech/DiNT)\
Release download of DiNT v2.0 available at: [UMN](https://comp.chem.umn.edu/dint/)\
Manual in PDF format: [Manual](https://comp.chem.umn.edu/dint/Dint_v2.0_manual.pdf)

## Authors
DiNT v1.0 (Jul 4 2013): Ahren W. Jasper, C. Melania Oana, and Donald G. Truhlar
DiNT v2.0 (Dec 16 2020): Ahren W. Jasper, Rui Ming Zhang, and Donald G. Truhlar

## Functionality
DiNT is a classical trajectory program for computing adiabatic and nonadiabatic chemistry. The code has a variety of options for preparing initial conditions, for performing final state analyses, and for semiclassical nonadiabatic propagation. The user may supply their own analytic potential energy surface, use one of the previously developed surfaces included in the distribution, or perform direct dynamics calculations. Recent applications include collisional energy transfer and spin-forbidden kinetics calculations.

## Preferred citations for this code
 1. A. W. Jasper, C. M. Oana, and D. G. Truhlar, DiNT, July 2013.
 2. A. W. Jasper and D. G. Truhlar, "Non-Born–Oppenheimer molecular dynamics for conical intersections, avoided crossings, and weak interactions," in Conical Intersections: Theory, Computation, and Experiment, edited by W. Domcke, D. R. Yarkony, and H. Koppel (World Scientific, 2011), pp. 375-412 (chp. 10) or Adv. Ser. Phys. Chem. 17, 375-412 (2011).
    
## Install
  ```console
  user@login:$ cd src/
  user@login:$ make
  ```

(Note: cmake option currently under development for v2.0, not guaranteed to work)
