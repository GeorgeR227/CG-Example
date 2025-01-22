using Polynomials
using LinearAlgebra
using Krylov
using Plots

vals = [2,7]
A = diagm(vals)
b = zeros(2)
f(x,y) = 0.5*[x,y]'*A*[x,y] - b'*[x,y]

xs = range(-1, 1, length = 100)
ys = range(-1, 1, length = 100)

z = f.(xs', ys)
plt = contour(xs, y, z; fill=true, color=:turbo, levels=25)

x0 = [-0.5, 0.1]

r0 = b - A*x0

r = r0
x = x0
maxiter = 10
res_xs = [x0[1]]
res_ys = [x0[2]]

for iter in 1:maxiter 
  alpha = dot(r, r)/dot(r, A*r)
  x .+= alpha .* r
  push!(res_xs, x[1])
  push!(res_ys, x[2])

  r .= b .- A * x
end
scatter!(res_xs, res_ys)

# scatter(vals, [0,0])

# function get_ith_Tpoly(i)
#   coeffs = zeros(i+1)
#   coeffs[end] = 1
#   return Polynomials.ChebyshevT(coeffs)
# end

# T = get_ith_Tpoly(2)

# xs = range(-0.5, 7.5, length=500)

# top_input = (7 .+ 2 .- 2 .* xs) ./ (7 .- 2)
# denom = Polynomials.evalpoly((7 + 2) / (7 - 2), T, false)

# plot!(xs, Polynomials.evalpoly.(top_input, T, false) ./ denom)
