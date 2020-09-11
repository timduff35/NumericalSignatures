restart
R=QQ[X, Y, Z,W]
x = random QQ
y = random QQ
z = random QQ
Rot = 1 /(1+x^2+y^2+z^2) * (matrix{{1+x*x-(y*y+z*z), 2*(x*y-z), 2*(x*z+y)},{2*(x*y+z), 1+y^2-(x*x+z*z), 2*(y*z-x)},{2*(x*z-y), 2*(y*z+x), 1 +z*z -(x*x+y*y)}})
assert(Rot*transpose Rot == id_(QQ^3) and det Rot == 1)
t = random(QQ^3,QQ^1)
C1 = id_(QQ^3) | matrix{{0},{0},{0}}
C2 = Rot | t

-- sphere
f = (X-3*W)^2 + (Y-3*W)^2 + (Z-3*W)^2 -W^2
contour1 = diff(W,f) * W
I = eliminate(ideal(f,contour1),W)
irrPrimes = decompose I
saturate(I,X^2+Y^2+Z^2) -- need to saturate "calibrating conic" X^2+Y^2+Z^2

world2 = transpose((Rot|t) * transpose vars R || matrix{{W}})
f2 = sub(f, world2)
center2 = gens ker C2
contour2 = (transpose center2 * jacobian f)_(0,0)
I2 = eliminate(ideal(f,contour2),W)

-- can we verify I ~ I2 up to SE(3)?
