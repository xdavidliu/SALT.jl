if !isdefined(:Salt)
	include("salt.jl")
end

#=========== parameters =================

N = [25, 20, 7];
Nc = 3;
M = [25, 20, 3];
eps = readdlm("eps3d.txt")
fprof = readdlm("fprof3d.txt")
Npml = [0, 0, 2];
LowerPML = 0;
BCPeriod = 0; nev = 1;
printnewton = 1;
h = [0.166666666666667, 0.173205080756888, .2];
wa = 1.7;
y = 2.0;
k = [0.0, 0.0, 0.0];
modenorm = 0.01;
ftol = 1.0e-7;
thresholdw_tol = 1.0e-7;


# =============== find first lasing mode passive resonance ==========


modeout = "pass3d";  
bl = [1, -1, 1];
wreal = 1.725
wimag = -0.0;  
Passive(N, Nc, M, Npml, LowerPML, bl, BCPeriod, 
h, wreal, wimag, k, modenorm, eps, fprof, modeout)



dD = 0.05;
Dmax = 0.15;
namesin = ["pass3d"];
namesout = ["after3d"];
Nm = 1;
Creeper(N, Nc, M, Npml, LowerPML, printnewton, 
h, wa, y, ftol, thresholdw_tol,
dD, Dmax, eps, fprof, namesin, namesout, Nm)



if !isinteractive()  # calling SlepcFinalize from an interactive session kills the session
	ccall((:SlepcFinalize, slepc), PetscErrorCode, ());
end