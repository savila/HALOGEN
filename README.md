=======
HALOGEN
=======

HALOGEN is a C program for creating realistic synthetic dark matter halo catalogs 
in a flash. 

It uses a statistical approach for generating the masses of the halos from a 
parameterised halo mass distribution, and placing them down in a spatial 
distribution with cheap 2LPT field as scaffolding (using a single-parameter
bias parameterisation).


Installation
------------
HALOGEN depends only on the libraries FFTW (v2?) and GSL. To set these and 
other specific options for compilation, edit the top part of the Makefile.defs
file in the top-level directory. Compilation has only been tested with mpicc/gcc.

Also in the Makefile.defs file are a number of global options for HALOGEN
which affect its output in various ways. Please set these as desired.

To make a version with in-built 2LPT capabilities, simple type ``make``.
To make a version which just performs the HALOGEN-specific components of
the calculation (which must be fed a 2LPT output), type ``make halogen``.
To make the fitting procedure, type ``make fit``.

The program will be compiled into the top-level directory.

To optionally install, save the executable somewhere on the system path, eg. 
in ``/usr/bin/``.


Usage
-----
There are 3 main programs that can be made, and each has a different
functionality.

The most comprehensive program is 2LPT-HALOGEN, which will perform a 
2LPT calculation, followed by a HALOGEN sampling on the resulting density
field. To run 2LPT-HALOGEN, type ``./2LPT-HALOGEN <2LPT.input>
<HALOGEN.input>``.

Example input files for each of these can be found in the examples/ 
subdirectory. Notably, the input files contain references to several
other files which must be provided. For the main HALOGEN input:

* alpha-file: A file containing a table of M and alpha, which encodes the 
  relative biasing of halos of different masses. For constant (with mass)
  bias, this should be one row long. This file can be produced with the
  ``fit`` program described below.
  
* hmf file: A file specifying the tabulated mass function with columns m,n(>m).
  This file can be generated in any way, but is most easily done with the 
  native support at hmf.icrar.org.
  
For the 2LPT input:

* glass file: A glass file

* Power Spectrum file: A table of logk,logP(k) for the desired cosmology. This may be
  produced in any fashion, but is simple with CAMB, CLASS, or at hmf.icrar.org.

If the standalone ``halogen`` program is used, type ``./halogen
<HALOGEN.input>``. In this case, an extra file needs to be included in the 
input specifications:

* 2LPT file: A gadget file, which is output of the 2LPT run. This could also
  be the output of a full NBODY run with gadget.

An important aspect of HALOGEN is a biasing relation parameterised as a
power-law of cell density. This relation may be defined in several halo mass
bins. For the purposes of reproducing simulation results, the optimal way
to determine this relation is by fitting directly to a realisation of the 
desired simulation for a specific setup. Note that this process is by far
the most time-consuming aspect of HALOGEN, but needs only to be done once for
a given cosmology, cell size and particle resolution. 

To run the fit, use ``??``.

  
Acknowledgments
---------------
If you find this code helpful in your research, please cite Avila, S. et
al. 2014 -- http://arxiv.org/abs/1412.5228.

Note that HALOGEN makes use of the pre-existing codes:
 - 2LPTic by Sebastian Pueblas and Roman Scoccimarro (http://arxiv.org/abs/astro-ph/9711187, http://cosmo.nyu.edu/roman/2LPT/)
 - CUTE  by David Alonso (http://arxiv.org/abs/1210.1833)

Authors
-------
Santiago Avila Perez: santiago.avila@uam.es
Steven Murray: steven.murray@uwa.edu.au 
