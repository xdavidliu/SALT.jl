# SALT

This is a solver under development for the SALT equations (steady-state
ab-initio lasing theory), solving them directly as a sparse nonlinear
system of equations, as described in:

* S. Esterhazy, D. Liu, M. Liertzer, A. Cerjan, L. Ge, K. G. Makris, A. D. Stone, J. M. Melenk, S. G. Johnson, and S. Rotter, “A scalable numerical approach for the steady-state ab-initio laser theory,” [Phys. Rev. A 90, 023816](http://journals.aps.org/pra/abstract/10.1103/PhysRevA.90.023816), August 2014.

It was originally written as a standalone C program, but we are
currently in the process of connecting it to a scripting front end in
the [Julia language](http://julialang.org/).
