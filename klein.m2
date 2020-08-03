restart
needs "main.m2"
setRandomSeed 0
dom = domain(3, 1)
Map = diffEqAffineSigMap dom
H = witnessHomotopy(dom, Map)
--Generic deg d EqAff Sig degree is 24d^2-48d (72 for d=3, 192 for d=4)
elapsedTime W = runMonodromy H

R = QQ[x,y,z]
--Should have Sig Deg = 18
f=x^3+y^3+z^3
--Should have Sig Deg = 178
f=x^3*y+y^3*z+z^3*x
Wf = witnessCollect(f, W,Verbose=>true)
