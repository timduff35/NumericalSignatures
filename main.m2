-*
Here we really just want the methods for witness collection and equality testing.

*-

needs "core.m2"
needs "mapDefinitions.m2"

-- core method: currently just does parameter homotopy from the witness data for a generic curve.
-- future methods should always return a WitnessData object
witnessCollect = method(Options=>{SlicePattern=>null,ImgTol=>1e-5, Verbose=>false, SliceParameters=>null, ChartParameters=>null})
witnessCollect (Matrix, WitnessData) := o -> (curveCoeffs, refW) -> (
    assert(numrows curveCoeffs == 1);
    numCoeffs := numcols curveCoeffs;
    H := homotopy refW;
    nSlcParams := numParameters gateSystem H - numParameters map H;
    nChartParameters := 4*dim domain H;
    chartParameters := if instance(o.ChartParameters, Nothing) then random(CC^1, CC^nChartParameters) else o.ChartParameters;
    sliceParameters := if instance(o.SliceParameters, Nothing) then random(CC^1,CC^nSlcParams) else o.SliceParameters;
    curveParameters := curveCoeffs | chartParameters | sliceParameters;
    trackWitness(H, refW, point curveParameters, ImgTol=>o.ImgTol, Verbose=>o.Verbose)
    )
witnessCollect (RingElement, WitnessData) := o -> (f, refW) -> (
    curveCoeffs := extractCoeffs f;
    witnessCollect(curveCoeffs, refW, o)
    )    

-*
Equality Test Input:
-- H: a WitnessHomotopy
-- Wref: WitnessData for the reference curve
-- testDomPt: representing a sample from the domain for the test curve
-- fixedParams: coefficients for the test curve and affine chart where the sample is defined (fixed by the homotopy)
Output:
  Boolean (are the curves equal up to symmetry?)...
    ... unless NearestWitnessPoint=>true: then we get the residual of the nearest reference witness point
        the latter is currently more expensive (no fast PointArray lookup)
        but still worth looking at in experiments where the runtime doesnt matter
*-

equalityTest = method(Options=>{
        NearestWitnessPoint => false, 
        Verbose=>false, 
        -- todo: these two options are possibly bit redundant. we should fix a convention for "bailing out" (currently we return null)
        SampleAttempts=>20, -- number of times we try to sample C
        TestAttempts=>5 -- number of test attempts w/ given sample from C
        })
equalityTest (RingElement, WitnessData) := o -> (f, refW) -> (
    C := extractCoeffs f;
    equalityTest(C, refW, o)
    )
equalityTest (Matrix, WitnessData) := o -> (C, refW) -> (
    H := homotopy refW;
    local tr;
    local result;
    attemptNum := 0;
    while (instance(result, Nothing) and attemptNum < o.TestAttempts) do (
        if o. Verbose then << "attempt number: " << toString(attemptNum) << endl;
        sampleResult := timing sampleCurveWParams(C, H,SampleAttempts=>o.SampleAttempts);
        trackTime := first sampleResult;
        (testDomPt, fixedParams) := last sampleResult;
        testCurveSliceParams := transpose sliceParams(point fixedParams, testDomPt, map H, slicePattern H);
        pRef := fixedParams |  sliceParams refW;
        pTest := fixedParams | testCurveSliceParams;
        tracking := timing trackWitness(H, point pTest, point pRef, {testDomPt}, Verbose=>o.Verbose);
        trackTime = trackTime + first tracking;
        imgPtRefWit := last tracking;
        print trackTime;
        theResult := if (length image imgPtRefWit == 0) then null else (
    	    imgPtRef := first image imgPtRefWit;
	    lookingUp := timing if (o.NearestWitnessPoint) then minPosition((points image refW)/(x->norm(matrix x-matrix imgPtRef))) else position(imgPtRef, image refW);
            << "^^^^ track time " << endl;
	    lookupTime := first lookingUp;
    	    matchedIndex := last lookingUp;
    	    testPositive := not instance(matchedIndex, Nothing);
    	    if testPositive then (
	        x := (image refW)_matchedIndex;
	        if o.NearestWitnessPoint then resid := norm(matrix x - matrix imgPtRef);
	        );
    	    if o.NearestWitnessPoint then resid else testPositive
    	    );
        tr = testResult(imgPtRefWit, theResult, trackTime, lookupTime);
        result = tr#"Result";
        attemptNum = attemptNum + 1;
        );
    tr
)


end--

--P2 diff Euclidean witness sets
restart
needs "main.m2"
setRandomSeed 2020
d=3
dom = domain(d,1)
Map = diffEuclideanSigMap dom
H = witnessHomotopy(dom,Map)

elapsedTime W = runMonodromy(H, NumberOfNodes=>3,Verbose=>true)

--teaser example
assert(d==3)
R = QQ[x,y,z]
f=homogenize(8*x^3 - (20*x)*y + 2*y^2 + 5*x - 10,z)
W0=witnessCollect(f,W,ImgTol=>1e-20,Verbose=>true)
length preimage W0

-- Klein quartic
assert(d==4)
R = QQ[x,y,z]
f=x^3*y+y^3*z+z^3*x
Wf = witnessCollect(f, W,Verbose=>true)

-- CHALLENGE (d=2)
R=QQ[x,y,z]
f = homogenize(20*x^2 + 35*x*y + 81*y^2 - 9*x - 8*y + 40,z)
Wf = witnessCollect(f, W,Verbose=>true)
evaluate(


restart
needs "main.m2"
setRandomSeed 0
d=3
dom = domain(d,4)
Map = jointEuclideanSigMap dom
H = witnessHomotopy(dom,Map)
elapsedTime W = runMonodromy(H, Verbose=>true)

R=QQ[x,y,z]

witnessCollect(f,W,Verbose=>true,ImgTol=>1e-20)
witnessCollect(f,W)

-- "bad

L0 = sub(matrix{{-100,60,3}},CC)
L1 = L0 + transpose foldVert for i from 1 to 1 list (matrix 0_(CC^2)||matrix{{random CC}})
L2 = L0 + transpose(matrix 0_(CC^2)||matrix{{random CC}})
IMGTOL=1e-10
W0 = witnessCollect(f, W,SliceParameters=>L0,ImgTol=>IMGTOL)
W1 = witnessCollect(f, W0,SliceParameters=>L1,ImgTol=>IMGTOL)
W2 = witnessCollect(f, W1,SliceParameters=>L2,ImgTol=>IMGTOL)
imSums = {W0,W1,W2}/(w->sum(matrix \ points image w))
r1=imSums#0-imSums#1
r2=imSums#0-imSums#2
realPoints points image W0
realPoints apply(preimage W0, x-> point evaluate(map H,point fixedParams W0,x))
SVD(r1||r2) -- trace test seems to pass
imSumIngredients = {W0,W1,W2}/(w->apply(preimage w, x -> evaluate(map H,point fixedParams w,x)))
imSumIngredients/first
imSums = imSumIngredients/sum
r1=imSums#0-imSums#1
r2=imSums#0-imSums#2
SVD(r1||r2) -- trace test seems to pass



-- Klein quartic
R = QQ[x,y,z]
f=x^3*y+y^3*z+z^3*x
setDefault(tStepMin=>1e-10)
--setDefault(
setDefault(maxCorrSteps=>1)
Wf = witnessCollect(f, W,Verbose=>true)
length preimage Wf -- 2712: todo, confirm this with trace test

-- trace test (on the previous random real curve)?
traceTest = method(Options=>{SVDTTolerance=>1e-6)
traceTest (WitnessData, WitnessHomotopy)
L0 = sliceParams W -- 28 = 4*7
L1 = L0 + transpose foldVert for i from 1 to 4 list (matrix 0_(CC^2)||matrix{{random CC}})
L2 = L0 + transpose(matrix 0_(CC^2)||matrix{{random CC}})
W0 = witnessCollect(f, W,SliceParameters=>L0)
W1 = witnessCollect(f, W0,SliceParameters=>L1)
W2 = witnessCollect(f, W1,SliceParameters=>L2)
imSums = {W0,W1,W2}/(w->sum(matrix \ points image w))
r1=imSums#0-imSums#1
r2=imSums#0-imSums#2
SVD(r1||r2) -- trace test seems to pass


outerIters = 10
innerIters = 10
setRandomSeed 0
k=infinity
eps = sub(10^(-k),RR)
C0 = sphere(binomial(d+2,2), 1);--random(RR^1,RR^(binomial(d+2,2)));
W0 = witnessCollect(C0, W, ImgTol=>5e-3, Verbose=>true);
C1 := randCoordChange(C0) + sphere(binomial(d+2,2),eps);    
timing equalityTest(C1,W0,ImgTol=>5e-3,NearestWitnessPoint=>(eps=!=0))

elapsedTime tests = for i from 1 to outerIters*innerIters list (
    if (-1+i%outerIters == 0) then (
	C0 = sphere(binomial(d+2,2), 1);--random(RR^1,RR^(binomial(d+2,2)));
    	W0 = witnessCollect(C0, W, ImgTol=>5e-3, Verbose=>true);
	);
    C1 := randCoordChange(C0,Noise=>eps);-- + sphere(binomial(d+2,2),eps);
    timing equalityTest(C1,W0,ImgTol=>5e-3
	,NearestWitnessPoint=>(eps=!=0)
	)
    );
iters = outerIters*innerIters
avgTime = sum(first \ tests)/iters
    nSucc = # select(last \ tests,status)
sum(select(tests / last,status) / getResult)/nSucc
netList(
 
    )

if areEual(eps, 0) then (

    nCorrect = # select(select(last \ tests,status), result)
    nCorrect / nSucc_RR) else (
    



netList sort apply(points image W,x->x)

x=symbol x
y=symbol y


R=CC[x,y,z]
f = random(3,R)
setDefault(tStepMin=>1e-11)
setDefault(ResidualTolerance=>1e-9)
--f=homogenize(8*x^3+6*x^2*y-7*x*y+2*y^2+5*x,z)
Wf = witnessCollect(f, W,Verbose=>true)
netList sort points image Wf
netList preimage Wf
pIm = homog \ points image Wf
netList(apply(pIm, x -> sort apply(pIm, y -> dP(x,y))))


pIm = homog \ points image W
netList(apply(pIm, x -> sort apply(pIm, y -> dP(x,y))))-- close!

x=symbol x
y=symbol y
R=RR[x,y,z]
f = random(3,R)

setDefault(tStepMin=>1e-11)
setDefault(ResidualTolerance=>1e-9)
--f=homogenize(8*x^3+6*x^2*y-7*x*y+2*y^2+5*x,z)
badSlice = matrix{{-1,2,1}}
goodSlice = matrix{{-5,1,10}}
worstSlice = matrix{{-2,0,1}}
Wf = witnessCollect(f, W,Verbose=>true, SliceParameters=>worstSlice)
netList sort points image Wf
netList preimage Wf
pIm = homog \ points image Wf
netList(apply(pIm, x -> sort apply(pIm, y -> dP(x,y))))

realPoints(
    (preimage Wf)/(x-> (
	    m := matrix x;
	    point((1/m_(0,2)) * m_{0,1})
	    )
    	)
    )
p := point fixedParams Wf
netList(
    {{"img","dom"}}|
    sort apply(preimage Wf, x -> {point evaluate(map H, p, x),dehomog x})
    
    

options witnessCollect

# preimage W
length points image W == 12 * (d^2-d)
netList sort preimage W
methods diff
D =diff(vars map H, gateMatrix map H);
Dgs = gateSystem(parameters map H, vars map H, transpose matrix{flatten entries D})
evD = (p, x) -> matrix(evaluate(Dgs,p,x),numrows D,numcols D)
netList sort apply(preimage W, x -> flatten entries evaluate(map H, point fixedParams W, x))
netList sort apply(preimage W, x -> (
	p := point fixedParams W;
	(flatten entries evaluate(map H, p, x), first SVD evD(p, x), x)
	    )
	) -- conditioning seems ok
netList sort apply(preimage W, x -> flatten entries evaluate(map H, point fixedParams W, point(10* matrix x)))
-- how close actually are these points?
netList apply(points image W, p -> sort apply(points image W, 
	q -> dP(
	    matrix{{1}}|matrix q,
	    matrix{{1}}|matrix p
	    )
	)
    )
-- this is a symptom of "high degree in low dimensions"

-- try random real curve: we observe the same issue for the map as previously defined
R=RR[x,y,z]
f = random(d,R) -- coefficients are sampled uniform [0,1]
Wreal = witnessCollect(f,W)
netList sort apply(preimage Wreal, x -> flatten entries evaluate(map H, point fixedParams Wreal, x))
--non-random real cuve?
Wreal = witnessCollect(y^3-z*y^2-x^3-x*z^2+z^3,W)
pIm = homog \ points image Wreal
netList(apply(pIm, x -> sort apply(pIm, y -> dP(x,y))))




length preimage Wreal
netList sort apply(preimage Wreal, x -> flatten entries evaluate(map H, point fixedParams Wreal, x))
netList sort apply(preimage Wreal, x -> (
	p := point fixedParams Wreal;
	(flatten entries evaluate(map H, p, x), first SVD evD(p, x), x)
	    )
	)

-- trace test (on the previous random real curve)?
L0 = sliceParams W
L1 = L0 + transpose(matrix 0_(CC^2)||matrix{{random CC}})
L2 = L0 + transpose(matrix 0_(CC^2)||matrix{{random CC}})
W0 = witnessCollect(f, W,SliceParameters=>L0)
W1 = witnessCollect(f, W0,SliceParameters=>L1)
W2 = witnessCollect(f, W1,SliceParameters=>L2)
imSums = {W0,W1,W2}/(w->sum(matrix \ points image w))
r1=imSums#0-imSums#1
r2=imSums#0-imSums#2
SVD(r1||r2) -- trace test seems to pass


# points image W == 12*(d^2-d)


-- k controls the noise
f = openOut "diffP2"
f << "noise, resid" << endl
tests = for k from -5 to -5 list new MutableList from {}
for k from -10 to -10 do (
    f << k;
    f << ", ";
    sig = 10^k;
    iters = 100;
    elapsedTime tests = for i from 1 to iters list (
	C0 := sphere(binomial(d+2,2), 1);--random(RR^1,RR^(binomial(d+2,2)));
	W0 := witnessCollect(C0, H, W);
	eps := 1e-3;
	C1 := randCoordChange C0 + sphere(binomial(d+2,2),eps);
	equalityTest(C1,W0,NearestWitnessPoint=>true)
    	);
    successes = select(tests,status);
    f << (1/(#successes))*sum(result \ successes) << endl;
    )
close f

-- w/ really big noise?
sig = 1
noisey = testCurveCoeffParams + sphere(10, sig)
-- we now need a way to sample the curve with those coefficients...
x3 = sampleCurve noisey
fixedParams = (noisey | transpose samples2ChartParams(x3))
equalityTest(H, W, samples2DomPt x3, fixedParams,NearestWitnessPoint=>true)


--joint Euclidean witness sets
-- gammify
restart
setRandomSeed 2020
needs "main.m2"
d=3
setRandomSeed 0
dom = domain(d,4)
Map = jointEuclideanSigMap dom
--SlicePattern = new HashTable from for i from 0 to 3 list {i} => 1 -- P1^6: cherry removed 
H = witnessHomotopy(dom,Map)
(p,x)=seedMonodromy H
elapsedTime W = runMonodromy(H,Verbose=>true)
d = degree domain H
D = dim domain H
q = point(gammaMat * transpose matrix p)
norm evaluate(H,p,x)
norm evaluate(H,q,x)

--H = witnessHomotopy(dom,Map,SlicePattern)
IMGTOL = 1e-5
elapsedTime W = runMonodromy(H,Verbose=>true,ImgTol=>IMGTOL)
# preimage W
length points image W
netList sort apply(preimage W, x -> flatten entries evaluate(map H, point fixedParams W, x))
outerIters = 10
innerIters = 10
setRandomSeed 0
elapsedTime tests = for i from 1 to outerIters*innerIters list (
    if (-1+i%outerIters == 0) then (
	C0 = sphere(binomial(d+2,2), 1);--random(RR^1,RR^(binomial(d+2,2)));
    	W0 = witnessCollect(C0, W, ImgTol=>5e-3, Verbose=>true);
	);
    C1 := randCoordChange C0;
    timing equalityTest(C1,W0,ImgTol=>5e-3)
    );
iters = outerIters*innerIters
avgTime = sum(first \ tests)/iters
nSucc = # select(last \ tests,status)
nCorrect = # select(select(last \ tests,status), result)
nCorrect / nSucc_RR -- only issues are due to failures in the underlying HC routine. I think this means the redefined map is invariant

outerIters = 10
innerIters = 10
setRandomSeed 0
iters = outerIters*innerIters
elapsedTime tests = for i from 1 to outerIters*innerIters list (
    if (-1+i%outerIters == 0) then (
	C0 = sphere(binomial(d+2,2), 1);--random(RR^1,RR^(binomial(d+2,2)));
    	W0 = witnessCollect(C0, H, W, ImgTol=>IMGTOL, Verbose=>true);
	);
    C1 := randCoordChange C0;
    timing equalityTest(C1,W0,ImgTol=>IMGTOL)
    );
avgTime = sum(first \ tests)/iters
nSucc = # select(last \ tests,status)
nCorrect = # select(select(last \ tests,status), result)
nCorrect / nSucc_RR


-- does witnessCollect work?
-- Fermat 
R=QQ[x,y,z]
f = x^d+y^d+z^d
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)
elapsedTime Wfermat = witnessCollect(f, H, W);
# points image Wfermat -- = 576 / 2, as expected

--P1xP1 diff Euclidean witness sets
restart
setRandomSeed 0
needs "main.m2"
d=10
setRandomSeed 0
dom = domain(d,1)
Map = diffEuclideanSigMap dom
SlicePattern = new HashTable from {{0}=>1}
H = witnessHomotopy(dom,Map,SlicePattern)
elapsedTime W = runMonodromy(H,Verbose=>true)


outerIters = 10
innerIters = 10
setRandomSeed 0
iters = innerIters*outerIters
elapsedTime tests = for i from 1 to outerIters*innerIters list (
    if (-1+i%outerIters == 0) then (
	C0 = sphere(binomial(d+2,2), 1);--random(RR^1,RR^(binomial(d+2,2)));
    	W0 = witnessCollect(C0, H, W, ImgTol=>5e-3, Verbose=>true);
	);
    C1 := randCoordChange C0;
    timing equalityTest(C1,W0,ImgTol=>5e-3)
    );
avgTime = sum(first \ tests)/iters
nSucc = # select(last \ tests,status)
nCorrect = # select(select(last \ tests,status), result)
nCorrect / nSucc_RR
-*
--P1xP1 diff Euclidean witness sets
d | deg
3 | 36
4 | 72
5 | 120
6 | 180
f = d -> 12*binomial(d,2)
f 6
*-
    
-*
--(P1)^5 joint Euclidean witness sets w/ cherry removed
d | deg(img) | timing
9 | 
8 | 
7 | 18816    | 1054 s
6 | 10080    | 384 s
5 | 4800     | 166 s
4 | 1920     | 47 s 
3 | 576      | 8 s
2 | 24       |
-- deg(d) = 8 d^2 (d^2-1)
*-




restart
-- projective invariants: generic cubic / quartic
needs "jointDemoProjective.m2"
d=3
setRandomSeed 0
(domGS, pts, paramValues, domPt) = setupDomain(d,9);
Map = projectiveSigMap(vars domGS, pts)
(F, p0, x0) = setupMonodromy(domGS, paramValues, domPt, Map, SlicePattern=> new HashTable from for i from 0 to 8 list {i} => 1) 
monodromySolve(F,p0,{x0},Verbose=>true)
Map = projectiveJointSigMap(domGS, pts)
--(F, p0, x0) = setupMonodromy(domGS, paramValues, domPt, Map)
out = 9
SLICES = remove(toList(0..9),out)
(F, p0, x0) = setupMonodromy(domGS, paramValues, domPt, Map, SlicePattern=> new HashTable from for i in SLICES list {i} => 1) 
setDefault(tStepMin=>1e-7)
setDefault(maxCorrSteps=>2)
elapsedTime monodromySolve(F,p0,{x0},Verbose=>true)

-*
--(P1)^5 Euclidean witness sets w/ disjoint edges removed
restart
needs "main.m2"
d=6
(domGS, pts, paramValues, domPt) = setupDomain(d,4);
Map = jointEuclideanSigMap(domGS, pts)
(F, p0, x0) = setupMonodromy(domGS, paramValues, domPt, Map, SlicePattern=> new HashTable from for i in {1,2,3,4} list {i} => 1) 
elapsedTime monodromySolve(F,p0,{x0},Verbose=>true)
d | deg(img) | timing
9 | 
8 | 
7 | 
6 | 13560    | 492 s
5 | 6304
4 | 2448
3 | 696      | 12 s
2 | 104       |
*-
d=symbol d
R=QQ[d]
rc = 11/6*d^4-40*d^3/3+164*d^2/3-817*d/6+141
2*rc

--Joint diff Euclidean witness sets
restart
needs "main.m2"
d=4
--setRandomSeed 0
right = 0
setDefault(tStepMin=>1e-8)
setDefault(maxCorrSteps=>2)
(domainGS, pts, paramValues, domPt) = setupDomain(d,2);
Map = jointDiffEuclideanSigMap(d, domainGS);
(F, p0, x0) = setupMonodromy(domainGS, paramValues, domPt, Map);
V = first monodromySolve(F,p0,{x0},Verbose=>true, NumberOfNodes=>3);
