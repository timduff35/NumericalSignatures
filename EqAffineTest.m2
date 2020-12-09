restart
needs "main.m2"
setRandomSeed 2
d=3
dom = domain(d,1)
Map = diffEqAffineSigMap dom
evaluate(Map, point random(CC^14,CC^1), point((random CC)*random(CC^3,CC^1)))


--generic signature has degree 72
H = witnessHomotopy(dom,Map)
elapsedTime W = runMonodromy(H, NumberOfNodes=>2,Verbose=>true)

--teaser example (should have deg 64 signature)
assert(d==3)
R = QQ[x,y,z]
f=homogenize(8*x^3 - (20*x)*y + 2*y^2 + 5*x - 10,z)
W0=witnessCollect(f,W,Verbose=>true)
length preimage W0

randCoordChange = method(Options=>{Group=>"E2"})
randCoordChange (Type, RingElement) := o -> (FF, f) -> randCoordChange(FF, first degree f, extractCoeffs f, o)
randCoordChange (Type, ZZ, Matrix) := o -> (FF, d, C) -> (
    M := randEqAff FF;
    transpose evaluate(
        coeffR d,
    	point sub((C_{0..binomial(d+2,2)-1}|matrix{flatten entries M}),CC)
	)
    )

equalityTest(randCoordChange(CC, f), W0)
