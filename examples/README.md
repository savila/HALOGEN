========
EXAMPLES
========
Included here are four examples of how to use HALOGEN. 

Each example is completely self-contained in its own directory.
To compile the relevant executables for these examples, edit
the top-level ``Makefile.defs`` according to your system
specifications (don't touch the ``DEFS`` at this stage), and
type ``make`` to make all the executables at once.

Each example is precisely the same where it overlaps with other
examples. Each is simulating a box called the GOLIAT simulation,
which is an L=1000, N=512 box with OmegaM = 0.27, 
OmegaLambda = 0.73, sigma_8 = 0.8 and H_0=70.0.

Running standalone 2LPT
-----------------------
HALOGEN includes a fully-working standalone 2LPT 
executable. This is included to be self-contained.
To compile the standalone 2LPT type ``make 2LPT``. 

An example input file for the ``2LPT`` executable is located 
in the ``ONLY_2LPT`` directory. To run it, use 

    $mpiexec -np 16 2LPT examples/ONLY_2LPT/GOLIAT_2LPT.input

from the top-level directory of HALOGEN. The number of mpi tasks
is specified after ``-np`` (16 in this example). But the specific
call for mpi will depend on the mpi compiler used.   

The input file for 2LPT includes settings for the boxsize and
number of particles, and the cosmology. In the case of this 
example, these are set to the GOLIAT simulation specs. 
It also includes paths to two extra files that are needed:

* *A glass file*: The one provided in the example is for a single
  particle species located at the centre of the grid.
* *A power spectrum file*: This is in the form of a table with two 
  columns, log10(k) and log10(P). Importantly, the power spectrum
  in this table includes the spectral tilt, n_s. Such a table may
  be simply downloaded using the inbuilt functionality at 
  hmf.icrar.org.

As set up, the output of the run will go into the ``output/``
subdirectory, and will contain a density field snapshot,
in GADGET format, at z=0.

Running standalone HALOGEN
--------------------------
HALOGEN can be built to run directly on a GADGET-format snapshot.
In this way, one can run many realisations of the HALOGEN-unique
part of the process on a single density field. It also serves to 
study how halogen parameters affect the output without changing
the density field. And it can also be used with other density 
fields different from 2LPT.

To compile the standalone HALOGEN executable, type ``make 
halogen``.

To run the example, located in the ``ONLY_HALOGEN`` directory,
make sure you first run the standalone 2LPT example above, as it is
the output from that run which will be used as the density field
for this example. 

Running the example is done simply with 

    ./halogen examples/ONLY_HALOGEN/GOLIAT_HALOGEN.input 

from the top-level directory of HALOGEN. The input for standalone 
HALOGEN, apart from the input specification file itself, consists of 
3 files:

* *Snapshot*: A GADGET-format snapshot of a density field. This can be
  in either GADGET 1 or 2 format. In the example, this is set to use
  the output of the previous 2LPT example, which is the GOLIAT box.
* *MassFunctionFile*: A tabulated cumulative mass function, with columns
  {M [M_sun/h], n(>M) [(h/Mpc)^3]}. The simplest way to get a mass 
  function is to use a fit from the literature. This can be easily done
  with the native functionality at http://hmf.icrar.org (make sure the cosmology
  is set appropriately). However, one could just as well use a function
  measured directly in a simulation.
* *alphaFile*: A table of mass bins, their corresponding alpha values and the velocity
  factor f_{vel}.
  This is the relation that defines linear bias and the velocity bias, respectively, in
  HALOGEN. The file can be produced by a fitting routine (see fitting example). The format
  of the file should be {alpha, M [M_sun/h], f_{vel}}, where M is the lower limit of
  each mass bin.
  Note that there is an example file at ``examples/FULL_HALOGEN/data/M-alpha.txt``, but if you ran ``./fit``
  you may also use the output of this ``examples/FITTING/output/GOLIAT/M-alpha.txt``. 

Besides this input data, the parameters defined in the HALOGEN.input file are
(values for this example in [square brackets]):

* *GadgetFormat*: [1]. Either 1 for Gadget 1 snapshot format, or 2 for Gadget 2
  format.
* *OutputFile* [``examples/ONLY_HALOGEN/output/GOLIAT.halos``]. The output file
  for the halo catalog.
* *NCellsLin*: [250]. The number of cells per side of the volume used in the
  HALOGEN method. Recommended setting is such that the number of particles
  per cell is about 2^3. In principle, the higher this number is the better,
  however, setting it too high could be problematic if warning messages saying 
  "MAXTRIALS reached" start to emerge.
* *recalc_frac*: [1.0]. A technical parameter of the method which controls how
  many times the probabilities of placing haloes in cells is updated.
  When the total relative error (over all cells) on the probailities
  exceeds *recalc_frac*, the recalculation is trigered.  
  Setting to 1.0 will correspond to only updating at mass-bin boundaries,
  which has been shown to be completely adequate. THIS PARAMETER MAY BE 
  REMOVED IN FUTURE VERSIONS.
* *nthreads*: [16]. The number of OPENMP threads to use in the placement
  routine. Note that in some machines, setting this too high can make it 
  work slowlier. 
* *rho_ref*: [crit]. Whether to take the halo overdensity definition as compared
  to critical density, or matter density (takes values {crit,matter}). 
* *Overdensity*: [200]. The overdensity criterion for the definition of a halo,
  R = ((3 M)/(4\pi rho_ref *Overdensity))^(1/3). This affects (together with 
  rho_ref) the exclusion radius of each halo. Halos will not be placed within 
  each-other's exclusion radius (see also the ``-DX_EXCLUSION`` defines in 
  ``Makefile.defs``).
* *Mmin*: [2e-4]. Sets either the minimum halo mass (in [M_sun/h]) or halo number 
  density (in [Mpc/h]^{-3}) of the final halo catalogue. 
  Which of these is set depends on the preprocessor define ``-DNDENS``.
* *Seed*: [-1]. The random seed of the halo placement algorithm. This can be set
  to ensure reproducible results, or left at -1, in which case the actual seed
  is calculated based on the current time. If one requires several
  HALOGEN-only simulations with different seeds, be careful that the
  calculations are spread out enough to receive different seeds.
* *Gad**: The various GADGET formatting parameters can be changed if one uses
  snapshots of slightly different format (e.g. change of units, 64 bits ID, double precision).
   These will be correct in the example if the HALOGEN version of 2LPT is used. 

The output of the HALOGEN-only run is a single halo catalog, located at
``output/GOLIAT.halos``, and is a table of halo positions (x,y,z) [Mpc/h], 
velocities (vx,vy,vz) [Mpc/h], masses [M_{sun}/h] and radii [Mpc/h].

Running Full HALOGEN
--------------------
The most common usage of HALOGEN will be to merge the ``2LPT`` and ``halogen``
components into one step, mitigating the expensive IO of writing the snapshot
to file. This is possible by using the ``2LPT-HALOGEN``
executable. Type ``make 2LPT-HALOGEN`` to compile it.  

To run the example, type

    $ mpiexec -np 16 2LPT-HALOGEN examples/FULL_HALOGEN/GOLIAT_2LPT.input examples/FULL_HALOGEN/GOLIAT_HALOGEN.input

The number of mpi tasks is specified after ``-np`` (16 in this example). 
But the specific call for mpi will depend on the mpi compiler used. 

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
``FULL_HALOGEN`` subdirectory).


Fitting
-------

This routine optimises the parameters which control the linear bias (alpha) and velocity bias
(f_vel). 

For running ``fit`` you should have run ``2LPT`` first or have a snapshot ready to input to this routine. 

To run the example in ``FITTING`` type: 

``./fit examples/FITTING/fit_GOLIAT.input``

The files needed by ``fit`` are the ones needed by the standalone ``HALOGEN`` plus a reference halo catalog:

* *NbodyFile*: Reference halo catalog from an N-Body simulation. ``fit`` will compute the 2PCF of this catalog
  and fit *alpha* to match that bias. It will also compute *fvel* from the velocities here. The halos should
  be in descending order of mass, and have the format X[Mpc/h], Y[Mpc/h], Z[Mpc/h], Vx[km/s], Vy[km/s], Vz[km/s], M[M_{sun}/h]

The input parameters for ``fit`` are very similar to those for ``halogen``. 
In addition, the output of ``fit`` must be specified in:

* *OutputDir* All the output of ``fit`` will be writen here with default names. Make sure that this directory is created, and that 
  it is specific for this fit (otherwise, it may be overwritten). The most important file will be called "M-alpha.txt" with columns
  {M alpha fvel}.

All the variables in commom between ``fit`` and ``halogen`` should be kept the same
when running ``fit`` and subsequently ``halogen`` (or ``2LPT-HALOGEN``) with
the "M-alpha.txt" a of the former (several of them, if changed, would require
a re-fit specific to those parameters).
 
The other (numerical) specific variables for ``fit`` are properly explained in
the example input file.
