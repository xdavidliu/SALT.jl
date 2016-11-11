# SALT

This is a solver under development for the SALT equations (steady-state
ab-initio lasing theory), solving them directly as a sparse nonlinear
system of equations, as described in:

* S. Esterhazy, D. Liu, M. Liertzer, A. Cerjan, L. Ge, K. G. Makris, A. D. Stone, J. M. Melenk, S. G. Johnson, and S. Rotter, “A scalable numerical approach for the steady-state ab-initio laser theory,” [Phys. Rev. A 90, 023816](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.90.023816), August 2014.

It was originally written as a standalone C program, but we are 
currently in the process of connecting it to a scripting front end in
the [Julia language](http://julialang.org/).

## Example

Here is an example of a simple, 1d ring laser with uniform dielectric and
gain. For convenience, all of the necessary files have be prepared, and can 
be found in the `examples/` directory. To find the first lasing threshold of
the passive pole with frequency 6.3, simply do

    ./run_ring passive
    ./run_ring creeper

The first command, `passive`, uses a linear eigensolver to find the passive poles
near some given frequency, and then outputs them to MATLAB `.m` mode files. (The reason for the MATLAB format is that PETSc defaults to this.) These files
are in general completely compatible with Julia, and can be read using the
`include` command, as usual, provided that `%` characters (comments in MATLAB) are replaced with `#` characters (comments in Julia).

The second command, `creeper`, takes some given `.m` mode files and increments the pump strength to some given value, while solving the SALT equations for all of the intermediate pump strengths. In this case, it is incrementing the pump strength for passive (below-threshold) modes until it finds the threshold, and then it outputs the new mode files at the threshold pump strength. The `creeper` command is also be used to solve for lasing (above-threshold) modes, which we will discuss later.


Now we examine the main executable BASH script, `run_ring`, which contains all the parameters necessary to find the lasing
mode for the 1d ring laser. 

The first part of the file gives the geometrical specifications:

    GEO="-Nc 1 -Nx 100 -Ny 1 -Nz 1
    -Npmlx 0 -Npmly 0 -Npmlz 0 -LowerPML 1
    -hx 0.01 -hy 0.1 -hz 0.1
    -epsfile "eps_ring.txt" -epsIfile "epsI_ring.txt" -fproffile "f_ring.txt"
    -wa 6.3 -gamma 1.0 -manual_epspml 0 -output_epstilde 0";

First, `-Nc` gives the
number of components of the electric field. For TM-polarized modes, it is 1;
for TE it is 2, and for fully vectorial fields it is 3. Next, `-Nx`, `-Ny`,
and `-Nz` give the number of total pixels, including any perfectly matched layer (PML) pixels, which we will describe next. For 1d geometries, both `-Ny` and `-Nz` are set to 1, while for 2d geometries, `-Nz` is set to 1. The `-Npml` fields give the number of pixels of PML placed at the boundaries. Here, we have set all to zero, because in the 1d ring, loss is modelled with an imaginary part of the dielectric, which we will soon describe. The `-LowerPML` parameter determines whether to put a PML on both the upper and lower boundaries (if set to 1) or just the upper boundary (if set to 0). The `-hx`, `-hy`, and `-hz` parameters determine the width (in arbitrary length units) of a single pixel. Note that our code chooses the speed of light to be 1, so these `h` parameters also directly determine the frequency units for all the results. Additionally, note that we have set both `-Ny` and `-Nz` to 1, so the values of `-hy` and `-hz` can be chosen to be anything, and do not matter. 

Next, we look at the input files. First, `-epsfile` must be set to an existing
text file that has the values of real parts of the dielectric function at all N×Ny×Nz grid points. In our case, we simply have "1" repeated 100 times in `eps_ring.txt`, since our ring has a uniform dielectric with value unity. The next file is `-epsIfile`, which specifies the imaginary part of the dielectric. Here, we have a text file `epsI_ring.txt` which has "0.2" repeated 100 times, to simulate a uniform radiation loss. For cavities that use outgoing boundary conditions to model radiation loss more realistically, we would instead have '-Npml' to be nonzero, and the `-epsIfile` option can simply be omitted, which will result in the code taking the dielectric to be purely real. The last file option, `-fproffile`, gives the file for the gain profile. Here we again have all ones. This file, unlike `-epsIfile`, is *not* optional.

 
