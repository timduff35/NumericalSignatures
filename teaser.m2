--teaser example
restart
needs "main.m2"
setRandomSeed 0
dom = domain(3, 1)
Map = diffEuclideanSigMap dom
H = witnessHomotopy(dom, Map)
elapsedTime W = runMonodromy H

R = QQ[x,y,z]
f=homogenize(8*x^3 - (20*x)*y + 2*y^2 + 5*x - 10,z)
W0 = witnessCollect(f,W)

-- intersect our signature curve with the blue line
-- notice in the output that we can't distinguish two points in the image due to rounding error
L0 = sub(matrix{{-100,60,-3}},CC)
W0 = witnessCollect(f, W,SliceParameters=>L0,Verbose=>true)

-- take slices parallel to L0 and run the trace test
-- we compare the result using 47 vs 48 points
L1 = L0 + transpose foldVert for i from 1 to 1 list (matrix 0_(CC^2)||matrix{{1_CC}})
L2 = L0 + transpose(matrix 0_(CC^2)||matrix{{-1_CC}})
W1 = witnessCollect(f, W0,SliceParameters=>L1)
W2 = witnessCollect(f, W1,SliceParameters=>L2)
imSums = {W0,W1,W2}/(w->sum(matrix \ points image w))
r1=imSums#0-imSums#1
r2=imSums#0-imSums#2
-- iffy as to whether or not the actual trace varies linearly based on our approximations to the points in the image (it does not)
min first SVD(r1||r2) 

-- now try with 48
imSumIngredients = {W0,W1,W2}/(w->apply(preimage w, x -> evaluate(map H,point fixedParams w,x)))
imSumIngredients/first
imSums = imSumIngredients/sum
r1=imSums#0-imSums#1
r2=imSums#0-imSums#2
min first SVD(r1||r2) -- our measure of linearity decreases by 3 orders of magnitude