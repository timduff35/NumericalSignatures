-- Example 4.1
restart
needs "main.m2"
setRandomSeed 2020
dom = domain(4, 1);
Map = diffEuclideanSigMap dom;
H = witnessHomotopy(dom, Map);
elapsedTime W = runMonodromy H; -- typically 9 seconds

R = QQ[x,y,z]
f=x^4+y^4+z^4
elapsedTime Wf = witnessCollect(f, W,Verbose=>true)

--joint
dom = domain(4, 4)
Map = jointEuclideanSigMap dom
H = witnessHomotopy(dom, Map)
elapsedTime W = runMonodromy H
elapsedTime Wf = witnessCollect(f, W,Verbose=>true)
