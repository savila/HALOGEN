========
EXAMPLES
========

Included here are three examples of how to use HALOGEN.

The expected usual workflow will be to run a fit for a given
simulation (ie. do something like the included FITTING/ example),
and then run the full ``2LPT-HALOGEN``, which merges the 2LPT and 
HALOGEN steps, skipping the expensive intervening I/O. 

Each example is completely self-contained in its own directory,
and includes its own readme.

Just remember that the M-alpha here is for j1Gp_N512_z0.00_WN_l4


Running 2LPT
------------
You can compile a standalone 2LPTic executable from within
HALOGEN. This is included to be self-contained. 

To compile the executable type ``make 2LPT``. Included here
is an example input file for the ``2LPT`` executable, located
in the ONLY_2LPT directory. To run it, use 
``./2LPT examples/ONLY_2LPT/2LPT.input`` from the top-level
directory of HALOGEN. 
