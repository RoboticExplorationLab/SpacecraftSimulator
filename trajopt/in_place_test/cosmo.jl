using LinearAlgebra



function proj(x)
    # i think its scalar first

    s = x[1]
    v = x[2:end]

    if norm(v)<=-s
        error("first case of projection")
    elseif norm(v)<=s
        return [s;v]
    else
        return .5*(1 + s/norm(v))*[norm(v);v]
    end
end


# input P q A b

ρ = .1

σ = 1e-6

α = 1.6

n = 3
m = 4

P = diagm(ones(n))
q = zeros(n)

A = diagm(ones(n))
b = zeros(m)
x = randn(n)
y = randn(n)
ν = randn(n)
s = randn(m)

max_iters = 10

LS_A = factorize([(P + σ*I)  A';
                   A         -(1/ρ)*I])

for k = 1:max_iters




    s_tilde = s - (1/ρ)*(ν + y)
