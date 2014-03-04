Nhalf = 50;
onevec = ones(Nhalf);
eps = [2.25*onevec, onevec];
fprof = [onevec, 0*onevec];
Npml = [20, 0, 0];
LowerPML = false;
printnewton = 1;
h = 0.5/Nhalf;
wa = 15.0;
y = 3.0;
k = [0.0, 0.0, 0.0];
modenorm = 0.01;
ftol = 1.0e-7;
thresholdw_tol = 1.0e-7;
modeout = "pass14";  
bl = [1, -1, 1];
wreal = 14.7;
wimag = -1.07;
geo = SALT.Geometry(eps, h, Npml, fprof, wa, y, n_vectorcomp=1, lowerPML=false);