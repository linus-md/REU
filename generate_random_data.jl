using HomotopyContinuation
using LinearAlgebra
using DataFrames
using CSV

@var s, t, r # Circle centered at (s, t) with radius r
@var u[1:3], v[1:3] # Three points of tangency (u[1], v[1]), ..., (u[3], v[3])
@var a[1:6], b[1:6], c[1:6] # Three fixed conics with coefficients defined by a[1], ..., a[6]; ... ; c[1], ..., c[6]

f_1 = (u[1] - s)^2 + (v[1] - t)^2 - r # the point (u[1], v[1]) lies on the circle
f_2 = (u[2] - s)^2 + (v[2] - t)^2 - r # the point (u[2], v[2]) lies on the circle
f_3 = (u[3] - s)^2 + (v[3] - t)^2 - r # the point (u[3], v[3]) lies on the circle
f_4 = a[1]*u[1]^2 + a[2]*u[1]*v[1] + a[3]*v[1]^2 + a[4]*u[1] + a[5]*v[1] + a[6] # the point (u[1], v[1]) lies on the conic defined by coefficients a[1:6]
f_5 = b[1]*u[2]^2 + b[2]*u[2]*v[2] + b[3]*v[2]^2 + b[4]*u[2] + b[5]*v[2] + b[6] # the point (u[2], v[2]) lies on the conic defined by coefficients b[1:6]
f_6 = c[1]*u[3]^2 + c[2]*u[3]*v[3] + c[3]*v[3]^2 + c[4]*u[3] + c[5]*v[3] + c[6] # the point (u[3], v[3]) lies on the conic defined by coefficients c[1:6]
f_7 = det([differentiate(f_1, [u[1], v[1]]) differentiate(f_4, [u[1], v[1]])]) # Circle and conic are tangent at (u[1], v[1])
f_8 = det([differentiate(f_2, [u[2], v[2]]) differentiate(f_5, [u[2], v[2]])]) # Circle and conic are tangent at (u[2], v[2])
f_9 = det([differentiate(f_3, [u[3], v[3]]) differentiate(f_6, [u[3], v[3]])]) # Circle and conic are tangent at (u[3], v[3])

paramVec = collect(Iterators.flatten([a,b,c]))
F = System([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9], variables = [u[1], v[1], u[2], v[2], u[3], v[3], s, t, r], parameters = paramVec);
df = DataFrame(A=Vector[], B=Int[])

for i in range(0, 1, 19500)
    p = randn(18)
    S = solve(F, target_parameters = p, show_progress = false);
    push!(df, [p, length(real_solutions(S))])
end

CSV.write("data_random.csv", df)
