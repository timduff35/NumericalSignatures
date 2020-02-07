restart
needs "main.m2"
setRandomSeed 0
dom = domain(4, 1)
Map = diffEuclideanSigMap dom
H = witnessHomotopy(dom, Map)
elapsedTime W = runMonodromy H

R = QQ[x,y,z]
--f=x^4+y^4+z^4
--f=x^3*y+y^3*z+z^3*x
Wf = witnessCollect(f, W,Verbose=>true)
