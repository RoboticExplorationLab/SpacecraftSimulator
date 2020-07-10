using LinearAlgebra



mutable struct tstruct
    a :: Float64
    B :: Array{Float64,2}
    C :: Array{Float64,2}
end



function runthetrap()

TT = tstruct(4.2,randn(3,3),randn(5,1000000))

a = 4.2
B = randn(3,3)
C = randn(5,1000000)

t1 = time()
for i = 1:1000000
    # a = randn()
    # B[rem(i,3)+1,rem(i,3)+1] = randn()
    # C[:,i] = randn(5)

    TT.a = randn()
    TT.B[rem(i,3)+1,rem(i,3)+1] = randn()
    TT.C[:,i] = randn(5)
end
t2 = time()
@show t2-t1

end

# runthetrap()

A = [2 1 1; 1 2 0; 1 0 2.0]

C = cholesky(A)


function chol_solve(L,b)
    """Solve (L*L')*x = b using forward and backwards substituion"""

    n = length(b)
    z = zeros(n)



    # L * (L'*x) = b
    # L * z      = b

    # forward substituion
    for i = 1:n
        α = b[i]

        for j = 1:(i-1)
            α = α - L[i,j]*z[j]
        end

        z[i] = α/L[i,i]
    end

    # back substition
    x = zeros(n)
    U = transpose(L)
    for i = n:-1:1
        α = z[i]

        for j = (i+1):n
            α = α - U[i,j]*x[j]
        end

        x[i] = α/U[i,i]
    end


    return z ,x
end
