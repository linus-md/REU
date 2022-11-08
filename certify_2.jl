using HomotopyContinuation
using LinearAlgebra
using DataFrames
using CSV

# The system F
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

F = System([f_1, f_2, f_3, f_4, f_5, f_6, f_7, f_8, f_9], 
                variables = [u; v; s; t; r],
                parameters = [a; b; c])

# Loading parameters
df = CSV.read("certify.csv", DataFrame)
params = map(1:69) do i
    p = Vector(df[i, 1:18])
end
numbers = map(1:69) do i
    no = df[i, 19]
end



# Solving
pt = randn(ComplexF64, 18)
St = solve(F, target_parameters = pt, show_progress = false)
S = solve(F, St, start_parameters = pt, target_parameters = params)

# Certification
C = map(S) do s
    res = s[1]
    p = s[2]
    certify(F, solutions(res), target_parameters = p)
end

# Compute number of strict complex solutions
c = certificates.(C)
number_of_strict_complex = map(c) do ct
    no = map(ct) do cert
        if !is_certified(cert)
            return false
        end
        tt = certified_solution_interval(cert)
        any(ttt -> !Base.in(0, imag(HomotopyContinuation.IComplexF64(ttt))), tt)
    end
    count(no)
end

# Compare
k = findall(nreal_certified.(C) .!= numbers)
t = findall((184 .- number_of_strict_complex) .!= numbers)