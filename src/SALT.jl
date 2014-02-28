module SALT

export Salt, Passive, Creeper

#=========== libraries ============ #
const saltlib = Pkg.dir("SALT", "deps", "saltlib");
const slepc = joinpath(ENV["SLEPC_DIR"], get(ENV, "PETSC_ARCH", ""), "lib", "libslepc");

typealias PetscErrorCode Cint;
const PETSC_COMM_WORLD = unsafe_load(cglobal((:PETSC_COMM_WORLD,slepc), Ptr{Void}));
ccall((:SlepcInitialize, slepc), PetscErrorCode,
              (Ptr{Cint}, Ptr{Ptr{Ptr{Uint8}}}, Ptr{Uint8}, Ptr{Uint8}),
              C_NULL, C_NULL, C_NULL, C_NULL);


#=========== master Salt function from C code ============ #

Salt(N, Nc, M, Npml, LowerPML, bl, BCPeriod, printnewton, 
h, wa, y, wreal, wimag, k, modenorm, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, modeout, namesin, namesout, Nm) = 
ccall((:Salt, saltlib), Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, 
Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cdouble, Cdouble,
Cint, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Uint8}, 
Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}}, Cint, 
Cint), int32(N), int32(M), h, int32(Npml), int32(Nc), int32(LowerPML), 
eps, fprof, wa, y, int32(BCPeriod), int32(bl), k, wreal, wimag, modenorm, 
int32(nev), modeout, dD, Dmax, thresholdw_tol, ftol, namesin, namesout, 
int32(printnewton), int32(Nm));


#============ sub-function for Passive resonance eigensolver ========== #

Passive(N, Nc, M, Npml, LowerPML, bl, BCPeriod, 
h, wreal, wimag, k, modenorm, eps, fprof, modeout) =
 
Salt(N, Nc, M, Npml, LowerPML, bl, BCPeriod, 0, 
h, 0.0, 0.0, wreal, wimag, k, modenorm, 1.0e-7, 1.0e-7,
0.0, 0.0, eps, fprof, modeout, [""], [""], 0);

#============ sub-function for Newton solver ======= #

Creeper(N, Nc, M, Npml, LowerPML, printnewton, 
h, wa, y, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, namesin, namesout, Nm) =

Salt(N, Nc, M, Npml, LowerPML, [0,0,0], 0, printnewton, 
h, wa, y, 0, 0, [0, 0, 0], 1.0, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, "", namesin, namesout, Nm);

include("plotSalt.jl")

end
