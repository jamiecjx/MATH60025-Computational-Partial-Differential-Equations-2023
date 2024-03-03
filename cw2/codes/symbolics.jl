using Symbolics
@variables h x(..) y(..)

xξ = (x(1,0) - x(-1,0))/(2h)
xη = (x(0,1) - x(0,-1))/(2h)
yξ = (y(1,0) - y(-1,0))/(2h)
yη = (y(0,1) - y(0,-1))/(2h)
xξξ =  (x(1,0) - 2x(0,0) + x(-1,0))/(h^2)
xξη =  (x(1,1) + x(-1,-1) - x(-1,1) - x(1,-1))/(4h^2)
xηη =  (x(0,1) - 2x(0,0) + x(0,-1))/(h^2)
yξξ =  (y(1,0) - 2y(0,0) + y(-1,0))/(h^2)
yξη =  (y(1,1) + y(-1,-1) - y(-1,1) - y(1,-1))/(4h^2)
yηη =  (y(0,1) - 2y(0,0) + y(0,-1))/(h^2)

α = xη^2 + yη^2
β = xξ*xη + yξ*yη
γ = xξ^2 + yξ^2


p1 = α * xξξ - 2β * xξη + γ * xηη




@variables u ux uxx uxy uy uyy uξ uξξ uξη uη uηη
@variables η ηx ηxx ηxy ηy ηyy 
@variables ξ ξx ξxx ξxy ξy ξyy 
@variables x xξ xξξ xξη xη xηη
@variables y yξ yξξ yξη yη yηη

@variables Jxy Jξη
Jξη = 1/ Jxy

ux = (yη * uξ - yξ * uη) / Jxy
uy = (yη * uξ - yξ * uη) / Jxy