needs "main.m2"
RL=symbol recursionLimit
MapAndRootCount = new HashTable from {
    "SE2" => (diffSEuclideanSigMap, d -> 6*d-d^2),
    "E2" => (diffEuclideanSigMap, d -> 12*d^2-12*d),
    "S2" => (diffSimilaritySigMap, d -> 9*d^2-13*d),
    "SA2" => (diffEqAffineSigMap, d -> 24*d^2-48*d),
    "A2" => (diffAffineSigMap, d -> 24*d^2-48*d),
    "PGL3" => (diffProjectiveSigMap, d -> 96*d^2-216*d)
    }

runWitnessExample = method(Options=>{Verbose=>true})
runWitnessExample (RingElement, HashTable) := o -> (f, CONFIG) -> (
    assert all(
        {Group, Degree, Seed, tStepMin, maxCorrSteps, RL, CorrectorTolerance}, 
        opt -> CONFIG#?opt
        );
    setDefault(tStepMin => CONFIG.tStepMin);
    setDefault(CorrectorTolerance => CONFIG.CorrectorTolerance);
    setDefault(maxCorrSteps => CONFIG.maxCorrSteps);
    recursionLimit := CONFIG.recursionLimit;
    d := CONFIG.Degree;
    assert(first degree f == d);
    G := CONFIG.Group;
    dom := domain(d, 1);
    (mapF, rcF) := MapAndRootCount#G;
    degIm := rcF d;
    Map := mapF dom;
    H := witnessHomotopy(dom, Map);
    elapsedTime W := runMonodromy(H,Verbose=>o.Verbose);
--    if o.Verbose then print length W.
    Wf := witnessCollect(f, W,Verbose=>o.Verbose);
    if o.Verbose then print Wf;
    (W, Wf)
    )
end
restart
needs "klein.m2"    

R = QQ[x,y,z]
f=x^3+y^3+z^3
--f=x^3*y+y^3*z+z^3*x

CONFIG = new HashTable from {
    Group => "A2",
    Degree => first degree f,
    Seed => 0,
    tStepMin => 1e-7,
    maxCorrSteps => 2,
    RL => 500, -- recursion limit
    CorrectorTolerance => 1e-6
    }

runWitnessExample(f, CONFIG)
(W,Wf) = oo
