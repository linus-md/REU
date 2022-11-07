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

p = append!(append!([0.01433106341596, 0.00768532229662, 0.00101592168808, 0.07448072207865, 0.02302883734631, -0.06521985689001, 0.01193897760052, 0.00377702785543, 0.00024606308255, 0.05128270773541, 0.00478610945491, 0.00263599776422, 0.01450050824971, 0.01393990808417, 0.00229893060745, 0.23894160029314, 0.05122052672597, 0.02150137716289]))
    #reverse([0.5965843419665127, -0.14037342740087336, -0.02041192362695606, 0.008257279526406407, 0.0024016615483861115, 0.0001737721803219947]), reverse([0.5965843418586586, -0.07399275257786551, 0.19629469624611, 0.0022925383444325925, -0.012173871420249097, 0.016146634896226787])), reverse([0.5965843418789551, 0.20701001428260404, -0.16136474402593104, 0.017957448375845633, -0.027997927138556127, 0.01090843408310055]))

function get_num(q)
    paramVec = collect(Iterators.flatten([a,b,c]))
    F = System(subs([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9], paramVec => q), variables = [u[1], v[1], u[2], v[2], u[3], v[3], s, t, r])
    S = solve(F, show_progress = false)
    l = length(real_solutions(S))
    return l
end

df = DataFrame(A=Vector[], B=Int[])

for i in range(0, 1, 500)
    d = randn(18)
    q = p + 10^-8 * d/norm(d)
    num = get_num(q)
    println(num)
    if num in [92, 94, 96, 98, 100, 102]
        push!(df, [q, num])
    end
end

CSV.write("data_100_2.csv", df)
