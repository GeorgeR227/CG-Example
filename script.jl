using Polynomials
using LinearAlgebra
using Krylov
using Plots

vals = [2,7]
A = diagm(vals)

scatter(vals, [0,0])

function get_ith_Tpoly(i)
  coeffs = zeros(i+1)
  coeffs[end] = 1
  return Polynomials.ChebyshevT(coeffs)
end

T = get_ith_Tpoly(2)

xs = range(-0.5, 7.5, length=500)

top_input = (7 .+ 2 .- 2 .* xs) ./ (7 .- 2)
denom = Polynomials.evalpoly((7 + 2) / (7 - 2), T, false)

plot!(xs, Polynomials.evalpoly.(top_input, T, false) ./ denom)
