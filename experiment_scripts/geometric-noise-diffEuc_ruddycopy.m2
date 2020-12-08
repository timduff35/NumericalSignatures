--teaser example
path = append(path, "../")
needs "main.m2"
setRandomSeed 0
d = 6
dom = domain(d, 4) -- diff (d,1) vs. joint (d,4)
Map = jointEuclideanSigMap dom -- diff (diffEuclideanSigMap) vs. joint (jointEuclideanSigMap)
H = witnessHomotopy(dom, Map)
-- phase 0 -- get witness set for signature of "generic curve"
NNODES=if d==2 then 3 else 2
elapsedTime W = runMonodromy(H, NumberOfNodes=>3)

-- R = QQ[x,y,z]
-- f=homogenize(8*x^3 - (20*x)*y + 2*y^2 + 5*x - 10,z)

-- phase 1: get witness set for Sig(f)
-- W0 = witnessCollect(f,W)
runExperiment = method(Options=>{FF=>CC, Verbose=>false, SubSampleSize=>null})
runExperiment (RingElement, WitnessData, RR) := o -> (f, W0, eps) -> (
    R := ring f;
    --1) get samples (x1,y1)...(xk,yk) on C (degree d)
    --samples := preimage W0;
    samples := apply(entries transpose sampleCurve(binomial(d+2,2)+1,f), l -> point{l});
    --2) "noise them" (X1,Y1)...(XK,YK)
    noiseySamples := samples/(s->matrix s + eps * sphere(3,1))/point;
    numSubSamples := if instance(o.SubSampleSize, Nothing) then #noiseySamples else o.SubSampleSize;
    randSubset := toList(0..numSubSamples-1);
    noiseySamples = apply(toList randSubset, i -> noiseySamples#i);
    --3) find "curve of best fit" C' from (X1,Y1)...(XK,YK)
    B := basis(d,R);
    GS := gateSystem(gateMatrix{getVarGates R}, transpose gateMatrix{gatePolynomial \ (sort flatten entries B)});
    M := matrix for s in noiseySamples list {evaluate(GS,s)};
    (S,U,Vt) := SVD M;
    coeffsC' := Vt^{numrows Vt-1};
    WC' := witnessCollect(coeffsC', W);
    -- 4) transform C' -> C''
    -- note: make CC vs RR an option
    coeffsC'' := randCoordChange(CC, d, coeffsC');
    -*
    5) run EQ(C, C'') (using transformed samples from part 3)
    EQ(C1,C2) 
      IN: 
       * samples on C1 which give a witness set for Sig(C1)
       * single sample from C2
      OUT: 
       * true/false
    *-
    ans := equalityTest(coeffsC'', W0, SampleAttempts=>50, TestAttempts=>50, Verbose=>o.Verbose);
    result := ans#"Result";
    if o.Verbose then (
        print result;
        if instance(result, Nothing) then print "test failed" else if result then print "curves are equivalent" else print "curves are not equivalent";
        );
    result
    )
end--
restart
load "geometric-noise-diffEuc_ruddycopy.m2"

setRandomSeed 0
R = CC[x,y,z]

randCurve = (d) -> (
    coeffs := sphere(binomial(d+2,2),1);
    monom := flatten apply(d+1,j -> apply(j+1, i -> x^i*y^(j-i)*z^(d-j)));
    f := sum(apply(binomial(d+2,2), i -> coeffs_i_0 * monom_i));
    f)
curves = apply(10, i -> randCurve(d))
witnesses = apply(10, i -> witnessCollect(curves_i,W))

tally flatten apply(10, i-> apply(10, j -> runExperiment(curves_i, witnesses_i, 1e-3)))
tally flatten apply(10, i-> apply(10, j -> runExperiment(curves_i, witnesses_i, 1e-4)))
tally flatten apply(10, i-> apply(10, j -> runExperiment(curves_i, witnesses_i, 1e-5)))
tally flatten apply(10, i-> apply(10, j -> runExperiment(curves_i, witnesses_i, 1e-6)))
tally flatten apply(10, i-> apply(10, j -> runExperiment(curves_i, witnesses_i, 1e-7)))

apply(entries transpose sampleCurve(11, curves_0), l -> point{l})
	
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

