-*
example: comparing near-quadrifolia
w/ differential signatures
*-

--offline phase 1: get witness set for generic curve
needs "main.m2"
dom = domain(6, 1)
Map = diffEuclideanSigMap dom
H = witnessHomotopy(dom, Map)
setRandomSeed 2020-- ensures monodromy result is reproducible
elapsedTime W = runMonodromy H

-- cast of characters
R = QQ[x,y]
quad=(x^2+y^2)^3-3*x^2+y^2
nearQuad = quad + 1/10000
t0 = 13213553/321749873941
R0 = 1/(t0^2+1)*matrix{{(t0^2-1),-2*t0},{2*t0, t0^2-1}}
T0 = random(QQ^1,QQ^2)
quadSE3 = sub(quad, vars R * R0 + T0)
nearQuadSE3 = sub(nearQuad, vars R * R0 + T0)


-- put curves in P^2
R=QQ[gens R|{z}]
f1=homogenize(sub(quad,R),z)
f2=homogenize(sub(quadSE3,R),z)
f3=homogenize(sub(nearQuadSE3,R),z)

-- offline phase 2: parameter homotopy to f1
elapsedTime W1=witnessCollect(f1, W)

 -- online phase: are f1 and f2 E2-equivalent?
tally apply(10, i -> result equalityTest(f2,W1))
equalityTest(f2,W1,NearestWitnessPoint=>true)
 -- with noise: are f1 and f3 E2-equivalent?
equalityTest(f3,W1)
equalityTest(f3,W1,NearestWitnessPoint=>true)

-- GB comparison demonstrates f1 and f3 are not erally equivalent
R=QQ[r_1,r_2,r_3,r_4,t_1,t_2][gens ring quad]
F1=sub(quad,R)
F3=sub(nearQuadSE3,R)
rot = matrix{{r_1,r_2},{r_3,r_4}}
xySUB = rot * matrix{{x},{y}} + matrix{{t_1},{t_2}}
F1F3equivEquations = last coefficients(F3 - sub(F1,transpose xySUB))
S=coefficientRing R
I=ideal(sub(F1F3equivEquations,S)) --+ ideal(rot*transpose rot - id_(S^2), (det rot)^2 - 1)
gens gb I
