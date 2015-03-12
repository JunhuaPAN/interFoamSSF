This is a multiphase solver based on the CFD framework OpenFOAM (version 2.0.x). It's an extension of the standard interFoam solver with the SSF formulation (see Raeini, Blunt, Bijeljic 2012 `[1]`) to reduce spurious currents.

`[1]` http://dx.doi.org/10.1016/j.jcp.2012.04.011

Known issues:
  * Dambreak test case is not low Capillary number

A discussion of the code can be found on [cfd-online](http://www.cfd-online.com/Forums/openfoam-programming-development/69673-parasitic-currents-2.html).

For a version that runs on OF-2.1.x, see [clones](http://code.google.com/p/interfoamssf/source/clones).