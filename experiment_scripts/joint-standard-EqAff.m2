restart
setRandomSeed 0
path = append(path,"../")
needs "main.m2"
needs "defaults.m2"

RESULTSFILE = openOut("../experiment_data/joint-standard-EqAff-degree-"|toString(d))

dom = domain(d,6);
Group = "SA"
Map = jointEqAffSigMap dom;
H = witnessHomotopy(dom,Map);

monodromyResult = timing runMonodromy(H, 
    NumberOfNodes=>NNODES, NumberOfEdges=>NEDGES, Verbose=>VERBOSE, 
    ImgTol=>IMGTOL, NumberOfRepeats=>NREPEATS, Gammify=>GAMMIFY);

randCoordChange = method(Options=>{Group=>"E2"})
randCoordChange (Type, RingElement) := o -> (FF, f) -> randCoordChange(FF, first degree f, extractCoeffs f, o)
randCoordChange (Type, ZZ, Matrix) := o -> (FF, d, C) -> (
    M := randEqAff FF;
    transpose evaluate(
        coeffR d,
    	point sub((C_{0..binomial(d+2,2)-1}|matrix{flatten entries M}),CC)
	)
    )

NOISELEVELS = {0.0} | for k from -9 to -5 list 10.0^k

load "master-experiment-joint.m2";
close RESULTSFILE

end--
restart
for i from 3 to 6 do (
    d = i;
    << " on degree " << d << endl;
    elapsedTime load "joint-standard-EqAff.m2"
    )
