using Krylov
using LinearAlgebra
using Plots
using SparseArrays

vals = [0,7]
A = diagm(vals)
b = zeros(2)
f(x,y) = 0.5*[x,y]'*A*[x,y] - b'*[x,y]

xs = range(-1, 1, length = 100)
ys = range(-1, 1, length = 100)

z = f.(xs', ys)
plt = contour(xs, ys, z; fill=true, color=:turbo, levels=25)

x0 = [-0.5, 0]

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

b = rand(20)

A = diagm(vcat(5*ones(10), 100*ones(10)))
cg(A, b)

A = diagm(range(4.99, 5.01, length=20))
cg(A, b)

A = diagm(range(5, 100, length=20))
x, stats = cg(A, b)

x[1] += 100
x0 = x

x, stats = cg(A, b, x0)


A = diagm(range(99.99, 100.01, length=20))
cg(A, b)

A = diagm(range(5, 100, length=20))

Q,_ = qr(rand(20,20))
rotA = Q'*A*Q
eigvals(rotA)
cg(rotA, b)
