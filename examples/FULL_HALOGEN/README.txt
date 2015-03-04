====================
FULL HALOGEN EXAMPLE
====================

This folder contains a fully working example of running ``2LPT-HALOGEN``.
For simplicity, all necessary input files have been provided, and the 
simulation corresponds to a boxsize of 1Gpc, with 512^3 particles, in a 
flat Universe with OmegaM = 0.27, H0 = 70.0, and sigma_8 = 0.8 at 
redshift 0 (these parameters can be seen in the file 2LPT.input).

Running the Example
-------------------
To run the example, copy ``2LPT-HALOGEN`` to the directory in which
this readme resides, then simply type:

    $ ./2LPT-HALOGEN HALOGEN.input 2LPT.input

A file called "halos.out" should be created in the output directory which 
is a table of halo positions (x,y,z), velocities (vx,vy,vz), masses 
(Msol/h) and radii (Mpc/h). There should be 200,000 haloes in this file.

Understanding the Input
-----------------------
The most important parameters in the 2LPT.input file are the cosmological
parameters, along with the parameters which specify the number of particles
used, and the boxsize. In essence, the 2LPT.input file constructs the overall
properties of the mock simulation. Note that a change in these parameters
necessitates a change in the HMF and input power spectrum, as well as a 
re-fitted alpha-M relation (thus one cannot expect to have reliable results
in this example when modifying just these parameters). 

In this example, we have provided a tabulated power spectrum produced with
CAMB. However, it is possible to produce the spectrum on-the-fly with an
Eisenstein-Hu fit or Efstathiou fit, though these are not as accurate.

Several parameters have been commented in this example as they are not 
required at all.

In the HALOGEN.input file, note that the file paths are all specified as 
relative -- this is relative to the executable, and therefore should not
in general be trusted. Use full paths to avoid confusion. 

The MassFunctionFile setting provides a path to a tabulated mass function
which corresponds to the correct cosmology and redshift. This may be obtained
from an analytic fit [eg. from hmf.icrar.org] or directly from a simulation.

Other important parameters are 

* NCellsLin, which specifies the number of cells per side of the box (and 
  therefore the resolution scale of the construction), and which we 
  recommend to set so that Lbox/NCellsLin ~ 4.
* alphaFile, which specifies a file which tabulates the relation of alpha
  with mass. This file simultaneously fixes the number of mass bins used
  in the reconstruction. This file has been provided as fit with the 
  parameters set in the input files. Any parameters which affect the output
  correlation function necessitate a recalibration of the alpha-M relation.
* rho_ref and Overdensity specify the exclusion region of each halo, defined
  as R = ((3 M)/(4\pi rho_ref *Overdensity))^(1/3). Re-specification of these
  values is unlikely to affect the alpha-M relation.
* Mmin, which specifies either the minimum halo mass used, or the total
  halo number density of the final output. Which is used is switched by
  the Makefile definition NDENS. 
