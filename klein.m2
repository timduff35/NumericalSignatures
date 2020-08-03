restart
setRandomSeed 0 --Fermat cubic has Similarity SigDeg = 12, Klein Q has Similarity SigDeg = 74
setRandomSeed 2020 --Fermat cubic has correct Similarity SigDeg = 10, Klein Q has correct Similarity SigDeg = 73
needs "main.m2"
dom = domain(4, 1)
Map = diffAffineSigMap dom
H = witnessHomotopy(dom, Map)

--Generic deg d EqAff Sig degree is 24d^2-48d (72 for d=3, 192 for d=4)
--Generic deg d Aff Sig degree is also 24d^2-48d (72 for d=3, 192 for d=4)
--Generic deg d Similarity Sig degree is 9*d^2-14*d (39 for d=3, 88 for d=4)
elapsedTime W = runMonodromy(H,Verbose=>true)

R = QQ[x,y,z]
--Should have EqAff Sig Deg = 18
--Should have Affine Sig Deg = 2
--Should have Similarity Sig Deg = 10
f=x^3+y^3+z^3
--Should have EqAff Sig Deg = 178
--Should have Affine Sig Deg = 24
--Should have Similarity Sig Deg = 73
f=x^3*y+y^3*z+z^3*x
setDefault(tStepMin=>1e-7)
Wf = witnessCollect(f, W,Verbose=>true)
