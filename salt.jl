const saltlib = joinpath(pwd(), "saltlib");
const slepc = joinpath(ENV["SLEPC_DIR"], get(ENV, "PETSC_ARCH", ""), "lib", "libslepc");
typealias PetscErrorCode Cint;
const PETSC_COMM_WORLD = unsafe_load(cglobal((:PETSC_COMM_WORLD,slepc), Ptr{Void}));
ccall((:SlepcInitialize, slepc), PetscErrorCode,
              (Ptr{Cint}, Ptr{Ptr{Ptr{Uint8}}}, Ptr{Uint8}, Ptr{Uint8}),
              C_NULL, C_NULL, C_NULL, C_NULL);


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


N = [100, 1, 1];
Nc = 1;
M = [50, 1, 1];
eps = 2.25*ones( M[1]*M[2]*M[3], 1);
fprof = 1.0*ones( M[1]*M[2]*M[3], 1);
Npml = [20, 0, 0];
LowerPML = 0;
BCPeriod = -1; nev = 1;
printnewton = 1;
h = [.01, .1, .2];
wa = 15.0;
y = 3.0;
k = [0.0, 0.0, 0.0];
modenorm = 0.01;
ftol = 1.0e-7;
thresholdw_tol = 1.0e-7;
dD = 0.0;
Dmax = 0.0;
namesin = [""];
namesout = [""];
Nm = 0;




modeout = "pass14";  
bl = [-1, -1, 1];
wreal = 14.7
wimag = -1.07;  # find first lasing mode passive resonance
Salt(N, Nc, M, Npml, LowerPML, bl, BCPeriod, printnewton, 
h, wa, y, wreal, wimag, k, modenorm, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, modeout, namesin, namesout, Nm);

modeout = "pass16";  
bl = [1, -1, 1];
wreal = 16.7
wimag = -1.07;  # find first lasing mode passive resonance
Salt(N, Nc, M, Npml, LowerPML, bl, BCPeriod, printnewton, 
h, wa, y, wreal, wimag, k, modenorm, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, modeout, namesin, namesout, Nm);


dD = 0.05;
Dmax = 0.5;
namesin = ["pass14", "pass16"];
namesout = ["after14", "after16"];
Nm = 2;
# turn up pump to D = 0.5
Salt(N, Nc, M, Npml, LowerPML, bl, BCPeriod, printnewton, 
h, wa, y, wreal, wimag, k, modenorm, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, modeout, namesin, namesout, Nm);




ccall((:SlepcFinalize, slepc), PetscErrorCode, ());