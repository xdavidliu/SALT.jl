N = [25, 20, 7];
Nc = 3;
eps = readdlm( Pkg.dir("SALT", "examples", "eps3d.txt") ) # todo: add air pixels
eps = reshape(eps, N[1], N[2], N[3]);
fprof = readdlm( Pkg.dir("SALT", "examples", "fprof3d.txt") )
fprof = reshape(fprof, N[1], N[2], N[3]);

Npml = [0, 0, 2];
LowerPML = false;
bl = [1, -1, 1];
wreal = 1.725; wimag = 0.0;
printnewton = 1;
h = [0.166666666666667, 0.173205080756888, .2];
wa = 1.7;
y = 2.0;
k = [0.0, 0.0, 0.0];
modenorm = 0.01;
ftol = 1.0e-7;
thresholdw_tol = 1.0e-7;