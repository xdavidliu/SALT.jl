const saltlib = "/Users/daveliu/Documents/saltc/saltlib";
const slepc = joinpath(ENV["SLEPC_DIR"], get(ENV, "PETSC_ARCH", ""), "lib", "libslepc");
typealias PetscErrorCode Cint  # from petscsys.h
const PETSC_COMM_WORLD = unsafe_load(cglobal((:PETSC_COMM_WORLD,slepc), Ptr{Void}));
ccall((:SlepcInitialize, slepc), PetscErrorCode,
              (Ptr{Cint}, Ptr{Ptr{Ptr{Uint8}}}, Ptr{Uint8}, Ptr{Uint8}),
              C_NULL, C_NULL, C_NULL, C_NULL);


N = int32([100, 1, 1]);
Nc = int32(1);
M = int32([50, 1, 1]);
Npml = int32([20, 0, 0]);
LowerPML = int32(0);
bl = int32([1, -1, 1]);
BCPeriod = int32(-1); nev = int32(1);
printnewton = int32(1);
h = [.01, .1, .2];
wa = 15.0;
y = 3.0;
wreal = 15.7;
wimag = -1.07;
k = [0.0, 0.0, 0.0];
modenorm = 0.01;
ftol = 1.0e-7;
thresholdw_tol = 1.0e-7;
dD = 0.05;
Dmax = 0.5;

epsfile = "eps1d.txt";
fproffile = "fprof1d.txt";
modeout = "pass";
namesin = ["pass0", "pass1"];
namesout = ["after0", "after1"];
Nm = 2;

ccall((:Salt, saltlib), Void, (Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Cint, Cint, Ptr{Uint8}, Ptr{Uint8}, Cdouble, Cdouble,
Cint, Ptr{Cint}, Ptr{Cdouble}, Cdouble, Cdouble, Cdouble, Cint, Ptr{Uint8}, Cdouble, Cdouble, Cdouble, Cdouble, Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}}, Cint, Cint), N, M, h, Npml, Nc, LowerPML, epsfile, fproffile, wa, y,
	BCPeriod, bl, k, wreal, wimag, modenorm, nev, modeout,
	dD, Dmax, thresholdw_tol, ftol, namesin, namesout, printnewton, Nm);