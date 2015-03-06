========
EXAMPLES
========

Included here are four examples of how to use HALOGEN.

The expected usual workflow will be to run a fit for a given
simulation (ie. do something like the included FITTING/ example),
and then run the full ``2LPT-HALOGEN``, which merges the 2LPT and 
HALOGEN steps, skipping the expensive intervening I/O. 

Each example is completely self-contained in its own directory,
and includes its own readme.

Just remember that the M-alpha here is for j1Gp_N512_z0.00_WN_l4


Running pure 2LPT
-----------------
You can compile a standalone 2LPTic executable from within
HALOGEN. This is included to be self-contained. 

To compile the executable type ``make 2LPT``. Included here
is an example input file for the ``2LPT`` executable, located
in the ONLY_2LPT directory. To run it, use 
``./2LPT examples/ONLY_2LPT/2LPT.input`` from the top-level
directory of HALOGEN. 

TIMING??

The input file for 2LPT includes settings for the boxsize and
number of particles, and the cosmology. In the case of this 
example, these are set to L=1 Gpc/h, N=512^3, and a cosmology with
OmegaM = 0.27, OmegaLambda = 0.73, sigma_8 = 0.8 and H0=70.0. It also 
includes paths to two extra files that are needed:

* A glass file. The one provided in the example is for a single
  particle species located at the centre of the grid.
* A power spectrum file. This is in the form of a table with two 
  columns, log10(k) and log10(P). Importantly, the power spectrum
  in this table includes the spectral tilt, n_s. Such a table may
  be simply downloaded using the inbuilt functionality at 
  hmf.icrar.org.

As set up, the output of the run will go into the ``output/``
subdirectory, and will contain a single density field snapshot,
in GADGET format, at z=0.

Running standalone HALOGEN
--------------------------
HALOGEN can be built to run directly on a GADGET-format snapshot.
In this way, one can run many realisations of the HALOGEN-unique
part of the process on a single density field. To compile the
standalone HALOGEN executable, type ``make halogen``. 

To run the example, located in the ``ONLY_HALOGEN`` directory,
make sure you first run the pure 2LPT example above, as it is
the output from that run which will be used as the density field
for the example. 

Running the example is done simply by ``./halogen
examples/ONLY_HALOGEN/HALOGEN.input`` from the top-level directory
of HALOGEN. The input for standalone HALOGEN, apart from the input
specification file itself, consists of three files:

* Snapshot: A GADGET-format snapshot of a density field. This can be
  in either GADGET 1/2 or 3 format. In the example, this is set to use
  the output of the previous 2LPT example, and is a box of L=1000 and 
  N=512.
* MassFunctionFile: A tabulated cumulative mass function, with columns
  {M [M_sun/h], n(>M) [(h/Mpc)^3]}. The simplest way to get a mass 
  function is to use a fit from the literature. This can be easily done
  with the native functionality at hmf.icrar.org (make sure the cosmology
  is set appropriately). However, one could just as well use a function
  measured directly in a simulation.
* alphaFile: A table of mass bins and their corresponding alpha values.
  This is the relation that defines the biasing scheme in HALOGEN, and the
  file can be produced by a fitting routine (see next example). The format
  of the file should be {alpha, M [M_sun/h]}, where M is the lower limit of
  each mass bin. VELOCITIES??

Besides this input data, the parameters defined in the HALOGEN.input file are
(values for this example in [square brackets]):

* GadgetFormat: [1]. Either 1 for Gadget 1 or 2 format, or 2 for Gadget 3
  format.
* OutputFile [examples/ONLY_HALOGEN/output/GOLIAT.halos]. The output file
  for the halo catalog.
* NCellsLin: [250]. The number of cells per side of the volume used in the
  HALOGEN method. Recommended setting is such that the number of particles
  per cell is about 2^3. In principle, the higher this number is the better,
  however, setting too high will introduce too much shot noise.
* Beta: ?????
* vel_fact: ?????
* recalc_frac: [1.0]. A technical parameter of the method which sets how
  many times the probabilities of placing haloes in cells is updated. 
  Setting to 1.0 will correspond to only updating at mass-bin boundaries,
  which has been shown to be completely adequate. THIS PARAMETER MAY BE 
  REMOVED IN FUTURE VERSIONS.
* nthreads: [16]. The number of OPENMP threads to use in the placement
  routine.
* rho_ref: [crit]. Whether to take the halo overdensity definition as compared
  to critical density, or matter density (takes values {crit,matter}). 
* Overdensity: [200]. The overdensity criterion for the definition of a halo,
  R = ((3 M)/(4\pi rho_ref *Overdensity))^(1/3). This affects (together with 
  rho_ref) the exclusion radius of each halo. Halos will not be placed within 
  each-others exclusion radius (see also the -DX_EXCLUSION defines in Makefile.defs).
* Mmin: [2e-4]. Sets either the minimum halo mass or total number density of 
  the final halo catalogue. Which of these is set depends on the preprocessor
  define -DDENS.
* Seed: [-1]. The random seed of the halo placement algorithm. This can be set
  to ensure reproducible results, or left at -1, in which case the actual seed
  is calculated based on the current time. If one requires several
  halogen-only simulations with different seeds, be careful that the
  calculations are spread out enough to receive different seeds.
* Gad*: The various GADGET formatting parameters can be changed if one uses
  snapshots of slightly different format. These will be correct in the example
  if the HALOGEN version of 2LPT is used. 

The output of the HALOGEN-only run is a single halo catalog, located at
``output/GOLIAT.halos``, and is a table of halo positions (x,y,z), 
velocities (vx,vy,vz), masses (Msol/h) and radii (Mpc/h).

Running Full HALOGEN
--------------------
The most common usage of HALOGEN will be to merge the ``2LPT` and ``HALOGEN``
components into one step, mitigating the expensive IO of writing the snapshot
to file. This is possible by using the ``2LPT-HALOGEN``
executable. Compiling this executable is done simply by typing ``make``.

To run the example, type

    $ ./2LPT-HALOGEN 2LPT.input HALOGEN.input

Besides the two input specification files, ``2LPT-HALOGEN`` requires all of
the input files necessary for the standalone ``2LPT`` and ``halogen``
executables, except for the Snapshot, which is never written to file. These
are all contained in the ``FULL_HALOGEN/data`` directory to make things
self-contained.

The input files for the full example are very similar to their standalone
counterparts, modulo parameters specific to the IO of the Snapshot (ie. the
Snapshot, GadgetFormat and Gad* parameters in the HALOGEN.input file, and the
output file parameter in the 2LPT.input).

The output should also be the same as the ``halogen`` executable (but in the 
FULL_HALOGEN subdirectory).


Fitting
-------
TODO!

