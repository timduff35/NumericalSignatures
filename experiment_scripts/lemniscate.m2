path=append(path,"../")
load "main.m2"
setDefault(tStepMin=>1e-8)
R=CC[x,y,z]
f=homogenize(
    (x^2+y^2)^2-(x^2-y^2)
    ,z)
dom = domain(4, 4)
Map = jointEuclideanSigMap dom
H = witnessHomotopy(dom, Map)
elapsedTime W = runMonodromy(H, Verbose=>true)
WJOINT = witnessCollect(f, W)
dom = domain(4, 1)
Map = diffEuclideanSigMap dom
H = witnessHomotopy(dom, Map)
elapsedTime W = runMonodromy(H, Verbose=>true)
WDIFF = witnessCollect(f, W)
end--
restart
setRandomSeed 0
needs "lemniscate.m2"
nsamples = 100
-- "perfect" samples
elapsedTime samples = for i from 0 to nsamples-1 list (
    t := 2*pi*i/nsamples; 
    matrix{{sin t/(1+(cos t)^2), sin t * cos t / (1+(cos t)^2),1}}
    )--apply(entries transpose sampleCurve(nsamples, f), s -> matrix{s});--for i from 0 to nsamples-1 list (t := 2*pi*i/nsamples; matrix{{2*(sin t)^2 * cos t, 2* (cos t)^2 * sin t,1}})


-- "unnoticeably noisey" samples
eps = 1e-3
noiseySamples = apply(samples, s -> s + (eps*sphere(2,1)|matrix{{0}}))
dat = openOut "lemniscateDataNoisy"
apply(flatten \ entries \ noiseySamples, s -> dat << toString(s#0) | ", " | toString(s#1) | "\n")
close dat
B = basis(4,R);
GS = gateSystem(gateMatrix{getVarGates R}, transpose gateMatrix{gatePolynomial \ (sort flatten entries B)});
M = matrix for s in noiseySamples list {evaluate(GS,s)};
(S,U,Vt) = SVD M;
coeffsC' = Vt^{numrows Vt-1};
clean_(1e-1)(coeffsC'||extractCoeffs f||matrix{sort flatten entries B})
first SVD(coeffsC'||extractCoeffs f)
eq = sub((coeffsC' * transpose matrix{sort flatten entries B})_(0,0),{last gens R=>1})
toString eq

diffAvg=1/100*sum for i from 0 to 99 list result equalityTest(coeffsC', WDIFF, NearestWitnessPoint=>true);
jointAvg=1/100*sum for i from 0 to 99 list result equalityTest(coeffsC', WJOINT, NearestWitnessPoint=>true);
randJointAvg=1/100*sum for i from 0 to 99 list result equalityTest(random(RR^1,RR^15), WJOINT, NearestWitnessPoint=>true);
randDiffAvg=1/100*sum for i from 0 to 99 list result equalityTest(random(RR^1,RR^15), WDIFF, NearestWitnessPoint=>true);
diffAvg
randDiffAvg
jointAvg
randJointAvg



-- "noticeably noisey" samples
eps = 5e-2
noiseySamples = apply(samples, s -> s + (eps*sphere(2,1)|matrix{{0}}))
dat = openOut "lemniscateDataNoisier"
apply(flatten \ entries \ noiseySamples, s -> dat << toString(s#0) | ", " | toString(s#1) | "\n")
close dat
B = basis(4,R);
GS = gateSystem(gateMatrix{getVarGates R}, transpose gateMatrix{gatePolynomial \ (sort flatten entries B)});
M = matrix for s in noiseySamples list {evaluate(GS,s)};
(S,U,Vt) = SVD M;
coeffsC' = Vt^{numrows Vt-1};
clean_(1e-1)(coeffsC'||extractCoeffs f)
first SVD(coeffsC'||extractCoeffs f)

errorDepth = 0
diffAvg=1/100*sum for i from 0 to 99 list result equalityTest(coeffsC', WDIFF, NearestWitnessPoint=>true);
jointAvg=1/100*sum for i from 0 to 99 list result equalityTest(coeffsC', WJOINT, NearestWitnessPoint=>true);
randJointAvg=1/100*sum for i from 0 to 99 list result equalityTest(random(RR^1,RR^15), WJOINT, NearestWitnessPoint=>true);
randDiffAvg=1/100*sum for i from 0 to 99 list result equalityTest(random(RR^1,RR^15), WDIFF, NearestWitnessPoint=>true);
diffAvg
randDiffAvg
jointAvg
randJointAvg
