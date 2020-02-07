-*

Ultimately, anything we don't want a user to see should go here.

Organization:

--0) imports, global variables, & uncategorized service functions
--1) creating and manipulating GateSystems and related objects
--2) manipulating points, parameters, and matrices
--3) methods for random sampling
--4) main classes

classes to consider adding: WitnessHomotopy, WitnessData (whats actually used in the equality test), WitnessPreImageData (needed for anything else), SlicePattern (shared by all of the previous)
*- 

--0) imports, global variables, & uncategorized service functions

needsPackage "MonodromySolver"
--debug needsPackage "NumericalAlgebraicGeometry"
debug needsPackage "NAGtypes"

-- global symbols that let us create slicing parameters on the fly
sliceCounter = 0
SS = symbol SS

-- fold objects with a "||" method in a list (ie matrices along rows)
foldVert = method()
foldVert List := L -> if (#L == 0) then L else fold(L,(a,b)->a||b)

-- fold objects with a "|" method in a list (ie matrices along cols)
foldHor = method()
foldHor List := L -> if (#L == 0) then L else fold(L,(a,b)->a|b)

-- complex arccos
arccos = z -> -ii*log(z+sqrt(z^2-1))

-- inverse of d -> (d+2)(d+1)/2
nDenseCoeffs2Degree = numCoeffs -> floor((1/2)*(-3 + sqrt(9+8*(numCoeffs-1))));

-- hermitian distance on projective space
dP = (v,w) -> max(0, arccos abs((1/(norm(2,v)*norm(2,w)))*transpose v*conjugate w)_(0,0))

--we're doing on the order of 100 tests, so I think we should only report percentages up to place after the decimal
pct = (ndecimals, r) -> sub(round(ndecimals,sub(100*r,RR)),RR)

-- number of rows and columns of matrices
size GateMatrix := J -> (numrows J, numcols J)
size Matrix := J -> (numrows J, numcols J)

-- extract coefficients from a homogeneous form in 3 variables with a consistent ordering
extractCoeffs = f -> (
    R := ring f;
    assert(numgens R == 3);
    F := coefficientRing R;
    assert(instance(F, InexactFieldFamily) or char F == 0);
    d := first degree f;
    S := F[gens R]; -- enforces grevlex
    transpose sub(last coefficients(sub(f, S), Monomials => sort basis(d, S)), CC)
    )

-*
-- TEST
v = random(CC^3,CC^1)
w = (random CC) * v
assert areEqual(dP(v,w),0)
*-

-- dehomogenize a point [x:y:z] in the chart z=1
dehomog = x -> (
    m := matrix x;
    point((1/m_(0,2)) * m_{0,1})
    )


--1) creating and manipulating GateSystems and related objects

-- gates for small determinants
det2 = M -> M_(0,0)*M_(1,1)-M_(1,0)*M_(0,1)
det3 = M -> M_(0,0)*det2(M^{1,2}_{1,2})-M_(0,1)*det2(M^{1,2}_{0,2})+M_(0,2)*det2(M^{1,2}_{0,1})

-- "join" of two GateSystems (take all functions from both)
GateSystem || GateSystem := (P, Q) -> (
    allVars := unique( (flatten entries vars P) | (flatten entries vars Q) );
    allParams := unique( (flatten entries parameters P) | (flatten entries parameters Q) );
    gateSystem(
	gateMatrix{allParams},
	gateMatrix{allVars},
    	(gateMatrix P)||(gateMatrix Q)
	)
    )

-- sum of two GateSystems
GateSystem + GateSystem := (P, Q) -> (
    if (numFunctions P =!= numFunctions Q) then error "can only add GateSystems of the same shape";
    H := P || Q;
    gateSystem(parameters H, vars H, gateMatrix P + gateMatrix Q)
    )

-- take some functions from the GateSystem
GateSystem ^ List := (P, inds) -> gateSystem(parameters P, vars P, (gateMatrix P)^inds)

-- append some slices to a given GateSystem
sliceSystem = method(Options => {Affine => 0, Homog => 0})
sliceSystem (GateMatrix, GateSystem) := o -> (X, F) -> (
    if (o.Affine <= 0 and o.Homog <= 0) then error("you did not do the slice");
    F || foldVert for i from 0 to o.Affine + o.Homog - 1 list (
	m := if i < o.Affine then numcols X + 1 else numcols X;
	X' := if i < o.Affine then transpose(X | gateMatrix{{1_CC}}) else transpose X;
    	sliceParams := gateMatrix{for i from 0 to m-1 list (
	    	ret := SS_sliceCounter;
	    	sliceCounter = sliceCounter + 1;
	    	ret)};
    	gateSystem(sliceParams, vars F, sliceParams * X')
	)
    )
sliceSystem GateSystem := o -> F -> sliceSystem(vars F, F, o)

-- sub new values for variables
sub (GateMatrix, GateSystem) := (X, F) -> gateSystem(parameters F, X, sub(gateMatrix F, vars F, X))    

-- return: a GateSystem for plane curve with parametric coefficients (full support) and a monomial matrix
genericCurve = d -> (
    C := symbol C;
    monIndices := apply(subsets(d+2,2), s -> (s#0,s#1-s#0-1,d+1-s#1));
    coeffs := apply(monIndices, i -> declareVariable C_i);
    monMat := mons(d, monIndices);
    gateSystem(gateMatrix{coeffs}, vars monMat, gateMatrix{coeffs} * gateMatrix monMat)
    )

-- gateSystem whose variables are group elements, parameters are curve coefficients
-- gives the action of a 3x3 matrix on polynomials of degree d
coeffR = d -> (
--    scan({C, g,x,y,z},s->s=symbol s);
    nGrpCoords := 9; -- we always transform by a 3x3 matrix
    Rng := CC[
	apply(subsets(d+2,2), s -> C_(s#0,s#1-s#0-1,d+1-s#1))| -- plane curve coeffs
	toList(g_(1,1)..g_(3,3)) -- group matrix entries
	][x,y,z];
    paramVars := drop(gens coefficientRing Rng, -nGrpCoords);
    G := transpose genericMatrix(coefficientRing Rng, g_(1,1), 3,3);
    f := sum(
	apply(paramVars, g -> (
	    	(dx,dy,dz) := last baseName g;
	    	g*x^dx*y^dy*z^dz
	    	)
	    )
    	);
    newCoeffs := coefficients sub(f,transpose(G*(transpose vars Rng)));
    symbCoeffs := apply(reverse flatten entries last newCoeffs,
	c->sub(c, coefficientRing Rng));
    polySystem symbCoeffs
    )


-- gateSystem for evaluating monomial matrix
mons = method()
mons ZZ := d -> mons(d, apply(subsets(d+2,2), s -> (s#0,s#1-s#0-1,d+1-s#1)))
mons (ZZ, List) := (d, monIndices) -> (
    C := symbol C;
    xx := symbol xx;
    yy := symbol yy;
    zz := symbol zz;
    R := CC[xx,yy,zz];
    gateSystem(polySystem apply(monIndices, i -> xx^(i#0)*yy^(i#1)*zz^(i#2)))
    )

--2) manipulating points, parameters, and matrices

-- reshape domain GateSystem's variables in a form suitable for defining a Map
reshape GateMatrix := mapVars -> (
    n := numcols mapVars;
    assert(n%3==0);
    foldVert for i from 0 to sub(n/3,ZZ)-1 list mapVars_{3*i,3*i+1,3*i+2}
)

-- given "samples" xs (a matrix whose columns represent samples from the curve)
-- get an affine chart on which each sample is defined, and then squash chart parameters into a matrix
samples2ChartParams = xs -> foldVert for i from 0 to numcols xs-1 list (
	numericalKernel(transpose xs_{i} | matrix{{1}}, 1e-5) * random(CC^3, CC^1)
	)

-- given a degree d and "samples" xs (a matrix whose columns represent samples from the curve)
-- this function returns a matrix of parameters consisting of
-- 1) coefficients of a generic curves of degree d containing the respective sample from xs, 
-- 2) plus random charts where each sample is defined
samples2CurveParams = (d, xs) -> (
    fCoefSpace := numericalKernel(
	foldVert for i from 0 to numcols xs-1 list evaluate(mons d,transpose xs_{i}),
	1e-5);
    fCoefValues := fCoefSpace * random(CC^(numcols fCoefSpace), CC^1);
    chartParamValues := samples2ChartParams xs;
    fCoefValues || chartParamValues
    )

-- collapses xs (a matrix whose columns represent samples from some curve) into a point suitable for evaluation
samples2DomPt = xs -> point matrix(transpose xs,(numcols xs)*(numrows xs),1);

--given values for domain variables, gateSystem representing Map, and a slice pattern, produce a random slice through the image of the point under the map
sliceParams = method(Options=>{Homog=>false})
sliceParams (Point, Point, GateSystem, HashTable) := o -> (mapParamPt, domPt, Map, slicePattern) -> (
    paramMat := random(CC^0,CC^1);
    scanPairs(slicePattern, (grp, nSlc) -> (
	    imgPt := (evaluate(Map, mapParamPt, domPt))_grp;
	    if not o.Homog then imgPt = imgPt |id_(CC^1);
	    K := numericalKernel(imgPt,1e-5);
	    sliceParams := K*random(CC^(numcols K),CC^nSlc);
--	    for i from 0 to numcols sliceParams -1 do print norm(2,sliceParams_{i});
	    paramMat = paramMat || matrix(transpose sliceParams,nSlc*(numrows K),1);
	    )
	);
    paramMat
    )


--3) methods for random sampling

-- coefficients of curve after random change of coordinates
randCoordChange = method(Options=>{Group=>"E2"})
randCoordChange (Type, RingElement) := o -> (FF, f) -> randCoordChange(FF, extractCoefficients f, o)
randCoordChange (Type, Matrix) := o -> (FF, C) -> (
    M := randE2 FF;
    testCurveCoeffPS = coeffR d; -- evaluate this to get new coefficients after coordinate change
    transpose evaluate(
        coeffR d,
    	point sub((C_{0..binomial(d+2,2)-1}|matrix{flatten entries M}),CC)
	)
    )


-- sample points from plane curve
sampleCurve = method(Options=>{})
sampleCurve Thing := o -> C -> sampleCurve(1, C, o)
sampleCurve (ZZ, RingElement) := o -> (k, f) -> (
    coeffs := extractCoefficients f;
    sampleCurve(k, coeffs, o)
    )
sampleCurve (ZZ, Matrix) := o -> (k, coeffs) -> (
    nCoeffs := numcols coeffs;
    d := nDenseCoeffs2Degree nCoeffs;
    mon := mons d;
    domGS := gateSystem domain(d,1,{symbol xx, symbol yy, symbol zz});
    I := gateSystem(vars domGS, transpose vars domGS);
    GS := domGS || (sliceSystem(I,Affine=>1))^{3};
    xs := random(CC^3,CC^k);
    foldHor for i from 0 to k-1 list (
    	domPt := samples2DomPt xs_{i};
	startParamValues := samples2CurveParams(d, xs_{i});
    	mapParamPt := point startParamValues;
    	startSliceParams = sliceParams(point{{}}, domPt, I, new HashTable from {{0,1,2}=>1});
    	startParamValues = startParamValues || startSliceParams;
    	targetParamValues := transpose coeffs || random(CC^(numrows startParamValues - nCoeffs), CC^1);    
    	transpose matrix first trackHomotopy(
	    specialize(parametricSegmentHomotopy GS, 
		startParamValues||targetParamValues
	    	),
	    {transpose matrix domPt}
	    )
    	)
    )	  

/// TEST
restart
needs "core.m2"
coeffs=random(CC^1,CC^10)
pt = transpose sampleCurve coeffs
M=mons 3
assert areEqual(0, norm evaluate(gateSystem(vars M,coeffs * (gateMatrix M)),pt))
///


-- random E2 element acting on P^2...of determinant -1
randE2 = FF -> (
    p := random FF;
    p1 := random FF;
    p2 := random FF;
    d := 1/(p^2+1);
    i := random 2;
    Rt := matrix{
	{d*(p^2-1),-2*d*p,p1},
	{d*2*p,d*(p^2-1),p2},
	{0,0,1}
	};
    ref := matrix{{-1,0,0},{0,1,0},{0,0,1}};
    ref*Rt
    )

---rnorm helper
gaussCC = () -> (
    (u1,u2):=(random RR,random RR);
    sqrt(-2*log(u1))*cos(2*pi*u2)+ii*sqrt(-2*log(u1))*sin(2*pi*u2)
    )

-- random sample drawn from normal distriution N(mu, var^2) via box-mueller method
rNorm = (mu,var) -> mu+var*(realPart gaussCC())_CC

-- random sample from (n-1)-sphere with radius r
sphere = (n,r) -> (
    l:=apply(n,i->rNorm(0,1));
    matrix{r/norm(2,l)*l}
    )

--4) main classes

-*
classes to consider adding: WitnessHomotopy, WitnessData (whats actually used in the equality test), WitnessPreImageData (needed for anything else), SlicePattern (shared by all of the previous)
*-

-- return: a parametric gateSystem representing a product of generic plane curves
-- note: this type is a dump!
Domain = new Type of HashTable
domain = method(Options=>{})
domain (ZZ, ZZ) := o -> (d, numDomainFactors) -> domain(d, genericCurve d, numDomainFactors, o)
domain (ZZ, ZZ, ZZ, List) := o -> (d, numDomainFactors, L) -> domain(d, genericCurve d, numDomainFactors, L, o)
domain (ZZ, GateSystem, ZZ, List) := o -> (d, fGS, numDomainFactors, inSymbs) -> (
    mapVars := matrix{declareVariable \ flatten for i from 1 to numDomainFactors list for s in inSymbs list  s_i};
    pts := reshape mapVars;
    domGS := foldVert for i from 0 to numDomainFactors-1 list sliceSystem(sub(pts^{i},fGS),Affine=>1);
    new Domain from {
	"GateSystem" => domGS,
	"Degree" => d,
	"Dim" => numDomainFactors
	}
    )
domain (ZZ, GateSystem, ZZ) := o -> (d, fGS, numDomainFactors) -> domain(d, fGS, numDomainFactors, {symbol x, symbol y, symbol z})
domain (ZZ, ZZ, List) := o -> (d, numDomainFactors, L) -> domain(d, genericCurve d, numDomainFactors, L, o)

gateSystem Domain := D -> D#"GateSystem"
degree Domain := D -> D#"Degree"
dim Domain := D -> D#"Dim"
net Domain := D -> net("domain == product of " | dim D | " curves of degree " | degree D)

WitnessHomotopy = new Type of HashTable 
witnessHomotopy = method(Options=>{})
witnessHomotopy (Domain, GateSystem) := o -> (dom, Map) -> (
    domGS := gateSystem dom;
    dimIm := numVariables domGS - numFunctions domGS;
    SlicePattern := new HashTable from {toList(0..numFunctions Map-1)=>dimIm};
    witnessHomotopy(dom, Map, SlicePattern)
    )
witnessHomotopy (Domain, GateSystem, HashTable) := o -> (dom, Map, SlicePattern) -> (
    GS := gateSystem dom;
    scanPairs(SlicePattern, (grp, nSlc) -> GS = sliceSystem((transpose gateMatrix Map)_grp, GS, Affine => nSlc));
    new WitnessHomotopy from {
	"Domain" => dom,
	"Map" => Map,
	"MasterGS" => GS,
	"SlicePattern" => SlicePattern
	}
    )

domain WitnessHomotopy := o -> H -> (H#"Domain")
map WitnessHomotopy := o -> H -> H#"Map"
gateSystem WitnessHomotopy := H -> H#"MasterGS"
slicePattern = method()
slicePattern WitnessHomotopy := H -> H#"SlicePattern"
dimIm = method() -- this presumably equals the number of factors in the domGS
dimIm WitnessHomotopy := H -> sum values slicePattern H
net WitnessHomotopy := H -> net("Witness Homotopy for " | net(domain H))
evaluate (WitnessHomotopy, Point, Point) := (H, p, x) -> evaluate(gateSystem H, p, x)
vars WitnessHomotopy := H -> vars gateSystem H
parameters WitnessHomotopy := H -> parameters gateSystem H
gateMatrix WitnessHomotopy := H -> gateMatrix gateSystem H

sampleCurveWParams = method(Options=>{SampleAttempts=>1})
sampleCurveWParams (Matrix, WitnessHomotopy) := o -> (C, H) -> (
    goodSample := false;
    attempts := 0;
    k := dim domain H;
    while not goodSample and attempts < o.SampleAttempts  do (
	samples := sampleCurve(k, C);
	p0Mat := C | transpose samples2ChartParams samples;
    	domPt := samples2DomPt samples;
	-- unclear why, but we get points that map to the origin of the signature curve with some non-negligible frequency
	-- we use a conservative tolerance to avoid this!
	mapEval := evaluate(map H, point p0Mat, domPt);
	if (instance(DBG, ZZ) and DBG > 1) then << mapEval << " w norm " << norm mapEval << endl;
	nEval := norm mapEval;
	if (not (areEqual(nEval, 0, Tolerance=>1e-1) or nEval > 1e6)) then goodSample = true else attempts = attempts + 1
	);
    if not goodSample then << "WARNING: image of samples on signature curve lie near origin" << endl;
    (domPt, p0Mat)    
    )
///TEST
restart
needs "core.m2"
needs "mapDefinitions.m2"
dom = domain(3,1);
Map = diffEuclideanSigMap dom;
SlicePattern = new HashTable from {{0}=>1};
H = witnessHomotopy(dom,Map,SlicePattern);
coeffs=sphere(10,1)
DBG=2
sampleCurveWParams(coeffs,H,SampleAttempts=>20)
M=mons 3
assert areEqual(0, norm evaluate(gateSystem(vars M,coeffs * (gateMatrix M)),pt))
///

jacobianCheck = (G, p0, domPt) -> (
    << "checking Jacobian" << endl;
    J := diff(vars G,gateMatrix G);
    (M,N) := size J;
    JGS := gateSystem(parameters G, vars G, transpose matrix{flatten entries J});
    J0 := matrix(transpose evaluate(JGS,p0,domPt),M,N);
    first SVD J0
    )

seedMonodromy = method(Options=>{SampleAttempts=>1,JacobianCheck=>false,ResidualCheck=>false})
seedMonodromy WitnessHomotopy := o -> H -> (
    d := degree domain H;
    C := random(CC^1,CC^(binomial(d+2,2)));--sphere(binomial(d+2,2),2);--sphere(1,random(CC^1,CC^(binomial(d+2,2)));
    (domPt, p0Mat) := sampleCurveWParams(C, H, SampleAttempts=>o.SampleAttempts);
    p0Mat = p0Mat | transpose sliceParams(point p0Mat, domPt, map H, slicePattern H);
    p0 := point p0Mat;
    G := gateSystem H;
    if o.ResidualCheck then (
	<< "checking residual" << endl;
	if not areEqual(norm evaluate(H, p0, domPt), 0.0) then error "seed point evaluation nonzero";
	);
    if o.JacobianCheck then (
	S := jacobianCheck(G, p0, x0);
    	if areEqual(min S, 0.0) then error "jacobian not full rank";
	);
    (p0, domPt)
    )



-- it should know: 
-- a WitnessHomotopy
-- variable and parameter values for the GateSystems associated to that WitnessHomotopy
-- points in the preimage (in a List, in case of future tracking needs)
-- points in the image (in a PointArray, for fast lookup)
WitnessData = new Type of HashTable
witnessData = method(Options=>{ImgTol=>1e-5,UsePointArray=>true})
witnessData (PointArray, Point, WitnessHomotopy) := o -> (pts, pt, H) -> witnessData(points pts, matrix pt, H, o)
witnessData (List, Matrix, WitnessHomotopy) := o -> (L, p, H) -> (
    new WitnessData from {
	"Homotopy" => H,
	"Parameters" => p,
	"PreimageWitnesses" => L,
	"ImageTolerance" => o.ImgTol,
	"ImageWitnesses" => (
	    img := clusterSolutions(
		L/(x->point evaluate(map H, point p_{0..numParameters Map -1}, x)),
		Tolerance => o.ImgTol
		);
	    P := pointArray {};
	    for i from 0 to #img -1 do if not member(img#i, P) then appendPoints(P, {img#i});
	    P
	    )
    }
)	 
witnessData WitnessData := o -> W -> witnessData(preimage W, parameters W, homotopy W, o)

net WitnessData := W -> net("witness data w/ " | length image W | " image points (" | # preimage W | " preimage points)")

preimage WitnessData := W -> W#"PreimageWitnesses"
parameters WitnessData := W -> W#"Parameters"
image WitnessData := W -> W#"ImageWitnesses"
homotopy = method()
homotopy WitnessData := W -> W#"Homotopy"
sliceParams WitnessData := o -> W -> (parameters W)_{numParameters map homotopy W..numParameters gateSystem homotopy W-1}
fixedParams = method()
fixedParams WitnessData :=  W -> (parameters W)_{0..numParameters map homotopy W-1}

pCompose = method()
pCompose (MutableHashTable, MutableHashTable) := (H1, H2) -> (
    new MutableHashTable from apply(keys H1,k-> if H2#?(H1#k) then k=> H2#(H1#k))
    )

writePermutations = (L, filename) -> (
    perms := L/(P->P/(i->i+1)); -- increment letters by 1 for GAP
    file := openOut (currentFileDirectory | filename);
    for i from 0 to #perms-1 do file << "p" << i << ":= PermList(" << toString(new Array from perms#i) << ");" << endl;
    file << "G:=Group(";
    for i from 0 to #perms-2 do file << "p" << i << ", ";
    file << "p" << #perms-1 << ");";
    close file;
    )


--INPUT: WitnessHomotopy
--OUTPUT: WitnessData
runMonodromy = method(Options=>{
	ImgTol=>1e-5,Verbose=>false,NumberOfNodes=>2,NumberOfEdges=>4,Filename=>null,NumberOfRepeats=>10,TargetSolutionCount=>null,SampleAttempts=>10, 
	SelectEdgeAndDirection=>selectRandomEdgeAndDirection, Potential=>null, Gammify=>true}) 
runMonodromy WitnessHomotopy := o -> H -> (
    (p0, domPt) := seedMonodromy(H,SampleAttempts=>o.SampleAttempts);
    getPerms := not instance(o.Filename, Nothing);
    d := degree domain H;
    nCoeffParams := binomial(d+2,2);
    D := dim domain H; -- == codim domain H! ie, this gives both the number of charts and the number of slices
    nMapCoords := 1+numrows gateMatrix map H;
    SlicePattern := slicePattern H;
    -- by default, we randomize the endpoints of monodromy by applying a random torus action
    randomizer := if not o.Gammify then null else (
	p -> (
	    gammaMat := diagonalMatrix(
		toList(
    		    (nCoeffParams:(random CC)) | 
    		    (foldHor for i from 1 to D list 4:(random CC))
		    ) | 
    		flatten apply(pairs SlicePattern, (grp, nSlc) -> toList foldHor for i from 1 to nSlc list (1+#grp):(random CC))
    		);
	    gammaMat * matrix p
	    )
	);
    V := first monodromySolve(gateSystem H, p0, {domPt},
	Verbose=>o.Verbose, NumberOfNodes=>o.NumberOfNodes, NumberOfEdges=>o.NumberOfEdges, EdgesSaturated=>true,
	NumberOfRepeats=>o.NumberOfRepeats,TargetSolutionCount=>o.TargetSolutionCount,
	SelectEdgeAndDirection=>o.SelectEdgeAndDirection, Potential=>o.Potential,
	Randomizer=>randomizer);
    if getPerms then (
	print "writing monodromy group file";
    	G := V.Graph;
    	V1 := first G.Vertices;
    	V2 := last G.Vertices;
    	E1 := toList V1.Edges;
    	-- "petal loops" based at V1
    	e1 := first E1;
    	perms := apply(drop(E1,1),e->values pCompose(e1.Correspondence12,e.Correspondence21));
    	writePermutations(perms,o.Filename)
    	);
    VV := first V.Graph.Vertices;
    witnessData(VV.PartialSols, VV.BasePoint, H, ImgTol=>o.ImgTol)
    )

-- can be used to either get a new witness set for the same Map with different slice values
-- or to parameter-homotop to a new curve
trackWitness = method(Options=>{ImgTol=>1e-5, Verbose=>false})
trackWitness (WitnessHomotopy, WitnessData, Point) := o -> (H, Wstart, pTarg) -> trackWitness(H, point parameters Wstart, pTarg, preimage Wstart, o)
trackWitness (WitnessHomotopy, Point, Point, List) := o -> (H, pStart, pTarg, startSols) -> (
    targetPreimages := trackHomotopy(
	specialize(
	    parametricSegmentHomotopy gateSystem H,
	    transpose(matrix pStart | matrix pTarg)
	    ),
	startSols
	);
    classes := partition(x -> status x == Regular, targetPreimages);
    if o.Verbose and classes#?false then scan(classes#false, x -> << " path failed due to " << status x << " at timestep " << x.LastT << endl);
    if (instance(DBG,ZZ) and DBG>0 and classes#?false) then error("investigate failure");
    retList := if classes#?true then (
	if o.Verbose then << # classes#true << " successful paths out of " << #targetPreimages << endl;
	classes#true) else {};
    witnessData(retList, matrix pTarg, H, ImgTol=>o.ImgTol)
    )

TestResult = new Type of HashTable
testResult = method()
testResult (WitnessData, Thing, RR, Thing) := (W, result, trackTime, lookupTime) -> new TestResult from {
    "Status" => length image W =!= 0,
    "Result" => result,
    "TrackTime" => trackTime,
    "LookupTime" => lookupTime
    }
status TestResult := o -> tr -> tr#"Status"
result = method()
result TestResult := tr -> tr#"Result"
trackTime = tr -> tr#"TrackTime"
lookupTime = tr -> tr#"LookupTime"