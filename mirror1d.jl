if !isdefined(:Salt)
	include("salt.jl")
end

#=========== parameters =================

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


# =============== find first lasing mode passive resonance ==========


modeout = "pass14";  
bl = [-1, -1, 1];
wreal = 14.7
wimag = -1.07;  
Passive(N, Nc, M, Npml, LowerPML, bl, BCPeriod, 
h, wreal, wimag, k, modenorm, eps, fprof, modeout)


#================ find second lasing mode passive resonance ==========

modeout = "pass16";  
bl = [1, -1, 1];
wreal = 16.7
wimag = -1.07;  
Passive(N, Nc, M, Npml, LowerPML, bl, BCPeriod, 
h, wreal, wimag, k, modenorm, eps, fprof, modeout)


dD = 0.05;
Dmax = 0.5;
namesin = ["pass14", "pass16"];
namesout = ["after14", "after16"];
Nm = 2;
# turn up pump to D = 0.5
Creeper(N, Nc, M, Npml, LowerPML, printnewton, 
h, wa, y, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, namesin, namesout, Nm)



if !isinteractive()  # calling SlepcFinalize from an interactive session kills the session
	ccall((:SlepcFinalize, slepc), PetscErrorCode, ());
end