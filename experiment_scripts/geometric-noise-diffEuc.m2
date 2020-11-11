--teaser example
path = append(path, "../")
needs "main.m2"
setRandomSeed 0
dom = domain(3, 1)
Map = diffEuclideanSigMap dom
H = witnessHomotopy(dom, Map)
-- phase 0 -- get witness set for signature of "generic curve"
elapsedTime W = runMonodromy H

R = QQ[x,y,z]
f=homogenize(8*x^3 - (20*x)*y + 2*y^2 + 5*x - 10,z)
-- phase 1: get witness set for Sig(f)
W0 = witnessCollect(f,W)
runExperiment = method(Options=>{FF=>CC, Verbose=>false, SubSampleSize=>null})
runExperiment (RingElement, WitnessData, RR) := o -> (f, W0, eps) -> (
    R := ring f;
    --1) get samples (x1,y1)...(xk,yk) on C (degree d)
    samples := preimage W0;
    --2) "noise them" (X1,Y1)...(XK,YK)
    noiseySamples := samples/(s->matrix s + eps * sphere(3,1))/point;
    numSubSamples := if instance(o.SubSampleSize, Nothing) then #noiseySamples else o.SubSampleSize;
    randSubset := toList(0..numSubSamples-1);
    noiseySamples = apply(toList randSubset, i -> noiseySamples#i);
    --3) find "curve of best fit" C' from (X1,Y1)...(XK,YK)
    B := basis(3,R);
    GS := gateSystem(gateMatrix{getVarGates R}, transpose gateMatrix{gatePolynomial \ (sort flatten entries B)});
    M := matrix for s in noiseySamples list {evaluate(GS,s)};
    (S,U,Vt) := SVD M;
    coeffsC' := Vt^{numrows Vt-1};
    WC' := witnessCollect(coeffsC', W);
    -- 4) transform C' -> C''
    -- note: make CC vs RR an option
    coeffsC'' := randCoordChange(CC, 3, coeffsC');
    -*
    5) run EQ(C, C'') (using transformed samples from part 3)
    EQ(C1,C2) 
      IN: 
       * samples on C1 which give a witness set for Sig(C1)
       * single sample from C2
      OUT: 
       * true/false
    *-
    ans := equalityTest(coeffsC'', W0, Verbose=>o.Verbose);
    result := ans#"Result";
    if o.Verbose then (
        print result;
        if instance(result, Nothing) then print "test failed" else if result then print "curves are equivalent" else print "curves are not equivalent";
        );
    result
    )
end--
restart
eps = 1e-10
load "geometric-noise-diffEuc.m2"

witnesses = apply(10, i -> witnessCollect(f,W))
W0 = witnesses#0

-- how does sample size (used for regressing the "approximate curve") effect outcome near threshold 10^(-4) 
setRandomSeed 0
tally apply(10, i -> runExperiment(f, W0, 1e-4, SubSampleSize => 11))
setRandomSeed 0
tally apply(10, i -> runExperiment(f, W0, 1e-4, SubSampleSize => 20))
setRandomSeed 0
tally apply(10, i -> runExperiment(f, W0, 1e-4, SubSampleSize => 40))
setRandomSeed 0
tally apply(10, i -> runExperiment(f, W0, 1e-4, SubSampleSize => 47))
setRandomSeed 0
tally apply(10, i -> runExperiment(f, W0, 1e-4, SubSampleSize => 48))

-- try different witnesses sets
setRandomSeed 2020 
toList apply(1..7, i -> (
        eps := sub(10^(-i), RR);
        (eps, tally apply(length witnesses, j -> runExperiment(f, witnesses#j, sub(10^(-i), RR), Verbose=>true)))
                    )
        )

