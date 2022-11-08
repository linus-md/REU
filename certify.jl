using HomotopyContinuation
using LinearAlgebra
using DataFrames
using CSV

df = CSV.read("/Users/las/Documents/GitHub/REU/certify.csv", DataFrame)
CSV.write("certify.csv", df, append=false)

@var s, t, r 
@var u[1:3], v[1:3] 
@var a[1:6], b[1:6], c[1:6] 

f_1 = (u[1] - s)^2 + (v[1] - t)^2 - r 
f_2 = (u[2] - s)^2 + (v[2] - t)^2 - r 
f_3 = (u[3] - s)^2 + (v[3] - t)^2 - r 
f_4 = a[1]*u[1]^2 + a[2]*u[1]*v[1] + a[3]*v[1]^2 + a[4]*u[1] + a[5]*v[1] + a[6] 
f_5 = b[1]*u[2]^2 + b[2]*u[2]*v[2] + b[3]*v[2]^2 + b[4]*u[2] + b[5]*v[2] + b[6]
f_6 = c[1]*u[3]^2 + c[2]*u[3]*v[3] + c[3]*v[3]^2 + c[4]*u[3] + c[5]*v[3] + c[6]
f_7 = det([differentiate(f_1, [u[1], v[1]]) differentiate(f_4, [u[1], v[1]])]) 
f_8 = det([differentiate(f_2, [u[2], v[2]]) differentiate(f_5, [u[2], v[2]])]) 
f_9 = det([differentiate(f_3, [u[3], v[3]]) differentiate(f_6, [u[3], v[3]])])


for i in 1:69
    p = Vector(df[i, 1:18])
    no = df[i, 19]
    paramVec = collect(Iterators.flatten([a,b,c]))
    F = System(subs([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9], paramVec => p), variables = [u[1], v[1], u[2], v[2], u[3], v[3], s, t, r])
    S = solve(F, show_progress = false)
    println(length(real_solutions(S)))
    println(certify(F, S))
end
