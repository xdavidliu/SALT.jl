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

Then make sure the `CREEPER` variable is set as

    CREEPER="-in0 pass_ring0
    -out0 thresh_ring
    -dD 0.07 -Dmax -0.001 -output_deps 0
    -newtonf_tol 1e-8 -reuse_lu 0 -printnewton 1";

Then finally, do

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

The remaining parameters in the block above are `-wa` and `-gamma`, which are the atomic gain frequency and relaxation the polarization rate from the SALT theory, respectively. Finally, `-manualepspml` and `-output_epstilde` are more advanced features, for manually providing a PML function (instead of the one computed from the `-Npml` parameters), and for outputting the total dielectric function for debugging purposes.

The next block of parameters, invoked when calling `passive`, is:

    PASSIVE="-wreal 6.2 -wimag -0.6
    -bx 1 -by -1 -bz 1
    -kx 1e-15 -ky 0 -kz 0 -norm 0.01
    -passiveout pass_ring -nev 2 -st_type sinvert"; 

The first line gives the guess for the pole when the `passive` command is used to find poles. The `-nev` parameter at the bottom gives, roughly, the number of eigenpairs to find. In this case, we are looking for poles near 6.2-0.6i, and we have told the eigensolver to save approximately two eigenpairs. The next set of parameters, `-b`, give the mirror boundary conditions. These are only relevant when the `-LowerPML` parameter is set to zero, when we only want to simulate half of the computational cell in each direction. In that case, setting `-bx` to -1 would impose odd mirror boundary conditions at the lower end in the x-direction. Other possibilities include 1 for even, and 0 for Dirichlet (zero field). In our case, however, we have set the `-LowerPML` parameter to 1 since we are trying to simulate the full cell, so all boundaries default to 0 (Dirichlet). 

The next line contains `-k`, which are the Bloch wave vectors in all directions. These can be used, for example, to simulate lasers in periodic media such as photonic crystals. In our case, we have set the `-kx` wave vector to 1e-15, which is equivalent to having periodic boundary conditions in the x direction. The next parameter `-norm`, simply provides an overall scaling to the electric field eigenvector. Usually, this parameter need not be adjusted; it is only when numerical stability becomes an issue that any modifications to it are necessary. Finally, the last line contains the `-passiveout` parameter, which specifies the filename prefix for the output mode files when the `passive` command is used. For example, since we have set the value of this option to `pass_ring`, the output files when `passive` is used would be `pass_ring0.m`, `pass_ring1.m`, and so on.

The third block of parameters, used when calling `creeper`, is given by:

    CREEPER="-in0 pass_ring0
    -out0 thresh_ring
    -dD 0.07 -Dmax -0.001 -output_deps 0
    -newtonf_tol 1e-8 -reuse_lu 0 -printnewton 1";

Here, the `-in0` and `-out0` parameters give the filenames (without the extension), for the input modes and output modes. The input file needs to exist prior to calling `creeper`. The initial pump strength is automatically read from the input mode file. For a situation with multiple *lasing* modes, one can use `-in1` and `-out1` for a *second* lasing mode, and the same with a 2 suffix for a *third* lasing mode, and so on. The caveat is that all the input modes must be at the same pump strength, or else an error might occur.

The next line of parameters include `-dD`, which gives the pump strength increment for each solve of the SALT equation, and `-Dmax`, which gives the final desired pump strength when it is a positive number, and indicates that we would like to find a *threshold* if it is negative. The rest of parameters here are more advanced settings.

After running the two commands, `./run_ring passive` and `./run_ring creeper` listed above, we will end up with a file `./thres_ring_file.m`, which is the mode exactly at threshold. We take a closer look at this file. The first variable in this file is `psi`, and it contains all the information about the mode. This `psi` variable always has Nx×Ny×Nz×Nc×2+2 elements. The first 2 in this expression comes from the fact that the field is complex and has real and imaginary components. The second 2 comes from the additional two variables, which are the real part of the frequency, and the imaginary part (if negative and the mode is below threshold) or the mode amplitude (if positive and above threshold). The remaining variables in the file are self-explanatory, with the exception of `ifix`, which is simply the position used to normalize and fix the phase of the mode, and can usually be ignored. For our threshold file, we see that the very last element of `psi` is indeed zero, which has been set explicitly by the `creeper` routine when the threshold was found to within a certain tolerance.

As a final step, we now want to take this threshold mode, and increase the pump strength even further so that starts lasing. To do so, we modify the `CREEPER` block of the BASH script to

CREEPER="-in0 thresh_ring
-out0 lasing_ring
-dD 0.07 -Dmax 0.4 -output_deps 0
-newtonf_tol 1e-8 -reuse_lu 0 -printnewton 1";

Note that we have moved the `thresh_ring` value from `-out0` to `-in0`, because we now want to use the existing file as an input. From this file, we can see the threshold is at approxiately `D=0.200024003819003;`. Hence, we choose `-Dmax` to be about twice this, and then do

    `./run_ring creeper`

The lasing mode at pump strength D = 0.4 is then output into the file `lasing_ring_file.m`.
