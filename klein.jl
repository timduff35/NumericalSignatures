# load the package
using HomotopyContinuation
using Random
using StatsBase


@var x y z
# Klein's quartic has the maximum # of automorphisms (168) of any genus-3 curve
f = x^3*y+y^3*z+z^3*x

# 168 is also the degree of the affine signature curve associated to f
# which we now compute a (pseudo)witness set for
y1 = - differentiate(f,x) / differentiate(f,y)
y2 = differentiate(y1,x) + differentiate(y1,y) * y1
y3 = differentiate(y2,x) + differentiate(y2,y) * y1
y4 = differentiate(y3,x) + differentiate(y3,y) * y1
y5 = differentiate(y4,x) + differentiate(y4,y) * y1
y6 = differentiate(y5,x) + differentiate(y5,y) * y1
T41 = 3*y4*y2;
T42 = 5*y3*y3;
T4 = T41 - T42;
T51 = 9*y5*y2*y2;
T52 = 45*y4*y3*y2;
T53 = 40*y3*y3*y3;
T5 = T51-T52+T53;
T61 = 9*y6*y2*y2*y2;
T62 = 63*y5*y3*y2*y2;
T63 = 45*y4*y4*y2*y2;
T64 = 255*y4*y3*y3*y2;
T65 = 160*y3*y3*y3*y3;
T6 = T61-T62-T63+T64-T65;
T4sq = T4^2;
T4cu = T4sq*T4;
denom1 = 1/T4sq;
denom2 = 1/T4cu;

@var u[1:4]
F=System([z-u[4],f,u[1] * denom2 * (T5)^2 + u[2] * denom1 * T6 + u[3]],parameters=u[1:4])
S=fixed(F,compile=true)

N=10
nsols = Array{Int64}(undef,N)
for i = 1:N
    Random.seed!(i);
    nsols[i] = length(solutions(monodromy_solve(S)))
end

show(countmap(nsols))
# (133 => 1,168 => 6,154 => 2,2 => 1)
# Galois group is imprimitive, so it's natural that monodromy is sometimes stuck on lower root count
# can probably make calculation more robust w/ more parameters
