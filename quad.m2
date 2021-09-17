
restart
needsPackage "MonodromySolver"
pts = apply(4, i -> (
    xi = random CC; yi = random CC; zi = sqrt((xi^2+yi^2)^3/(4*xi^2*yi^2)); point{{xi,yi,zi}}))
VARIABLES=gateMatrix{toList vars(x_1,y_1,z_1,x_2,y_2,z_2,x_3,y_3,z_3,x_4,y_4,z_4)}
PARAMETERS=gateMatrix{{declareVariable p} | toList vars(a_(1,1)..a_(4,7)) | toList vars(b_(1,1)..b_(4,4))}
x0=fold(matrix \ pts, (a,b) -> a|b)
xs = apply({1,2,3,4},i->x_i/z_i)
ys = apply({1,2,3,4},i->y_i/z_i)
dists = apply(subsets(4,2), S -> (
	(i,j):=(first sort S,last sort S);
	(xs#i-xs#j)^2+(ys#i-ys#j)^2
	)
    )
D = gateSystem(VARIABLES, transpose gateMatrix{dists})
eqs = (
    apply({1,2,3,4}, i -> (x_i^2+y_i^2)^3-p*x_i^2*y_i^2*z_i^2) |
    apply({1,2,3,4}, i -> b_(i,1)*x_i+b_(i,2)*y_i+b_(i,3)*z_i+b_(i,4)) |
    apply({1,2,3,4}, i -> a_(i,7) + sum apply(dists, {1,2,3,4,5,6}, (d, j) -> a_(i,j)*d))
    )
sys=gateSystem(PARAMETERS,VARIABLES,transpose gateMatrix{eqs})
(p0,x0)=createSeedPair(sys,point x0)
equivalencer = x -> point evaluate(D,x)
setDefault(maxCorrSteps=>2)
setDefault(tStepMin=>1e-9)
setDefault(CorrectorTolerance=>1e-8)
setDefault(PointArrayTol=>1e-3)
VNP=monodromySolve(sys, p0, {x0}, 
    Verbose=>true, 
--    Equivalencer => equivalencer, 
    PointArrayTol=>1e-4, NumberOfNodes=>4, NumberOfEdges => 1, Randomizer => (x->matrix x))
# monodromyGroup (first VNP).Graph
monodromyGroup((first VNP).Graph,FileName=>"test.txt");

