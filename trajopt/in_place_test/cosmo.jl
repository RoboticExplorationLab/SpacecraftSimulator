using LinearAlgebra, Convex, SCS, SparseArrays, Infiltrator


struct QCQP_struct
    ρ::Float64
    σ::Float64
    α::Float64
    P::Array{Float64,2}
    q::Array{Float64,1}
    u_max::Float64
    A::SparseMatrixCSC{Float64,Int64}
    b::SparseVector{Float64,Int64}
    LS_A::SparseMatrixCSC{Float64,Int64}
    LS_sol::Array{Float64,1}
    x_k::Array{Float64,1}
    y_k::Array{Float64,1}
    s_k::Array{Float64,1}
    xtilde_kp1::Array{Float64,1}
    ν_kp1::Array{Float64,1}
    stilde_kp1::Array{Float64,1}
    s_kp1::Array{Float64,1}
end


function create_QCQP(nu,u_max)

    # solver settings
    ρ = .1
    σ = 1e-6
    α = 1.6

    # cost matrices
    P = zeros(nu,nu)
    q = zeros(nu)

    # A
    A = sparse([1:nu;nu+1],[1:nu;nu],[ones(nu);0])
    b = sparsevec([nu+1],[u_max])

    # lsA
    LS_A = sparse([zeros(nu,nu) A'; A -(1/ρ)*I])

    # temporary variables
    LS_sol = zeros(2*nu+1)
    x_k = zeros(nu)
    y_k = zeros(nu+1)
    s_k = zeros(nu+1)
    xtilde_kp1 = zeros(nu)
    ν_kp1 = zeros(nu+1)
    stilde_kp1 = zeros(nu+1)
    s_kp1 = zeros(nu+1)

    QCQP = QCQP_struct(copy(ρ)::Float64,
                       copy(σ)::Float64,
                       copy(α)::Float64,
                       copy(P)::Array{Float64,2},
                       copy(q)::Array{Float64,1},
                       copy(u_max)::Float64,
                       copy(A)::SparseMatrixCSC{Float64,Int64},
                       copy(b)::SparseVector{Float64,Int64},
                       copy(LS_A)::SparseMatrixCSC{Float64,Int64},
                       copy(LS_sol)::Array{Float64,1},
                       copy(x_k)::Array{Float64,1},
                       copy(y_k)::Array{Float64,1},
                       copy(s_k)::Array{Float64,1},
                       copy(xtilde_kp1)::Array{Float64,1},
                       copy(ν_kp1)::Array{Float64,1},
                       copy(stilde_kp1)::Array{Float64,1},
                       copy(s_kp1)::Array{Float64,1})

end

function create_convex(nu,P,q,u_max)
    x = Variable(nu)
    problem=minimize(.5*quadform(x,P) + dot(q,x),[norm(x,2)<= u_max])

    return x, problem
end

function cone_proj(x)
    # i think its scalar first


    v = x[1:(end-1)]
    s = x[end]

    nv = norm(x[1:end-1])
    if nv<=-x[end]
        # @infiltrate
        # error()
        # println("first option")
        x .= zeros(length(x))
    elseif nv>= abs(s)
        x .= .5*(1 + x[end]/nv)*[v;nv]
    end

    return x
end


function QCQP_solve!(QCQP::QCQP_struct,P,q)

    # update problem matrices
    QCQP.P .= P
    QCQP.q .= q

    # update the linear system matrix
    QCQP.LS_A[1:length(q),1:length(q)] .= P + QCQP.σ*I

    # this is the new allocation each time
    # fA = factorize(QCQP.LS_A)
    fA2 = qdldl(QCQP.LS_A)
    nu = length(q)

    for k = 1:30


        # solve linear system
        # QCQP.LS_sol .= vec(fA\[(-QCQP.q + QCQP.σ*QCQP.x_k);(QCQP.b-QCQP.s_k + (1/QCQP.ρ)*QCQP.y_k)])
        QCQP.LS_sol .= solve(fA2,[(-QCQP.q + QCQP.σ*QCQP.x_k);(QCQP.b-QCQP.s_k + (1/QCQP.ρ)*QCQP.y_k)])


        QCQP.xtilde_kp1 .= QCQP.LS_sol[1:nu]
        QCQP.ν_kp1      .= QCQP.LS_sol[nu+1:end]

        # update stilde
        QCQP.stilde_kp1 .= QCQP.s_k - (1/QCQP.ρ)*(QCQP.ν_kp1 + QCQP.y_k)

        # update x
        QCQP.x_k .= QCQP.α*QCQP.xtilde_kp1 + (1-QCQP.α)*QCQP.x_k

        # update s
        QCQP.s_kp1 .= cone_proj(QCQP.α*QCQP.stilde_kp1 + (1-QCQP.α)*QCQP.s_k + (1/QCQP.ρ)*QCQP.y_k)

        # update y
        QCQP.y_k .= QCQP.y_k + QCQP.ρ*(QCQP.α*QCQP.stilde_kp1 + (1-QCQP.α)*QCQP.s_k - QCQP.s_kp1)

        # reset s
        QCQP.s_k .= QCQP.s_kp1

        if rem(k,5)==0
            if norm(QCQP.A*QCQP.x_k + QCQP.s_k - QCQP.b)<1e-6
                break
            end
        end

    end
end


function compare()

    nu = 3
    u_max = 3.67

    P = randn(nu,nu); P = P'*P;
    q = randn(nu)

    QCQP = create_QCQP(nu,u_max)

    # convex part
    x = Variable(nu)
    problem=minimize(.5*quadform(x,P) + dot(q,x),[norm(x,2)<= u_max])

    @time solve!(problem, SCS.Optimizer,verbose = false)

    @show x.value

    # my part of it
    @time QCQP_solve!(QCQP::QCQP_struct,P,q)

    @show QCQP.x_k

end


# function QC_solve(P,q,t_max)
#
#     # solver settings
#     ρ = .1
#     σ = 1e-6
#     α = 1.6
#     n = length(q)
#
#     A = [Diagonal(ones(n));zeros(n)']
#     A = sparse([1:n;n+1],[1:n;n],[ones(n);0])
#     b = sparsevec([n+1],[t_max])
#
#
#     LS_A = factorize(sparse([(P + σ*I)  A';
#                        A         -(1/ρ)*I]))
#     # LS_A =  inv([(P + σ*I)  A';
#     #                    A         -(1/ρ)*I])
#     LS_sol = zeros(n+m)
#     x_k = zeros(n)
#     y_k = zeros(m)
#     s_k = zeros(m)
#     xtilde_kp1 = zeros(n)
#     ν_kp1 = zeros(m)
#     for k = 1:max_iters
#
#
#         # solve linear system
#         LS_sol .= LS_A\[(-q + σ*x_k);(b-s_k + (1/ρ)*y_k)]
#         # LS_sol .= LS_A*[(-q + σ*x_k);(b-s_k + (1/ρ)*y_k)]
#         xtilde_kp1 .= LS_sol[1:n]
#         ν_kp1      .= LS_sol[n+1:end]
#
#         # update stilde
#         stilde_kp1 = s_k - (1/ρ)*(ν_kp1 + y_k)
#
#         # update x
#         x_k .= α*xtilde_kp1 + (1-α)*x_k
#
#         # update s
#         s_kp1 = cone_proj(α*stilde_kp1 + (1-α)*s_k + (1/ρ)*y_k)
#
#         # update y
#         y_k .= y_k + ρ*(α*stilde_kp1 + (1-α)*s_k - s_kp1)
#
#         # reset everything
#         # x_k .= x_kp1
#         # y_k .= y_kp1
#         s_k .= s_kp1
#
#         if rem(k,5)==0
#             err[k] = norm(x_k - x_true)
#             if err[k]<1e-6
#                 break
#             end
#         end
#
#     end


# function solveit()
# # input P q A b
#
# ρ = .1
#
# σ = 1e-6
#
# α = 1.6
#
# n = 3
# # m = 4
#
# P = randn(n,n);P = P'*P
# q = randn(n)
#
#
# t_max = 3.0
#
# A = [diagm(ones(n));zeros(n)']
# b = [zeros(n);t_max]
#
# m = size(A,1)
#
# # x_k = randn(n)
# # y_k = randn(m)
# # s_k = randn(m)
# x_k = zeros(n)
# y_k = zeros(m)
# s_k = zeros(m)
#
# max_iters = 30
#
#
#
# # LS_A = [(P + σ*I)  A';
# #          A         -(1/ρ)*I]
#
# # convex.jl part of it
# x = Variable(n)
# problem=minimize(.5*quadform(x,P) + dot(q,x),[norm(x,2)<= t_max])
#
# solve!(problem, SCS.Optimizer)
#
# x_true = copy(x.value)
#
# err = zeros(max_iters)
#
# t1 = time()
# LS_A = factorize(sparse([(P + σ*I)  A';
#                    A         -(1/ρ)*I]))
# # LS_A =  inv([(P + σ*I)  A';
# #                    A         -(1/ρ)*I])
# LS_sol = zeros(n+m)
# xtilde_kp1 = zeros(n)
# ν_kp1 = zeros(m)
# for k = 1:max_iters
#
#
#     # solve linear system
#     LS_sol .= LS_A\[(-q + σ*x_k);(b-s_k + (1/ρ)*y_k)]
#     # LS_sol .= LS_A*[(-q + σ*x_k);(b-s_k + (1/ρ)*y_k)]
#     xtilde_kp1 .= LS_sol[1:n]
#     ν_kp1      .= LS_sol[n+1:end]
#
#     # update stilde
#     stilde_kp1 = s_k - (1/ρ)*(ν_kp1 + y_k)
#
#     # update x
#     x_k .= α*xtilde_kp1 + (1-α)*x_k
#
#     # update s
#     s_kp1 = cone_proj(α*stilde_kp1 + (1-α)*s_k + (1/ρ)*y_k)
#
#     # update y
#     y_k .= y_k + ρ*(α*stilde_kp1 + (1-α)*s_k - s_kp1)
#
#     # reset everything
#     # x_k .= x_kp1
#     # y_k .= y_kp1
#     s_k .= s_kp1
#
#     if rem(k,5)==0
#         err[k] = norm(x_k - x_true)
#         if err[k]<1e-6
#             break
#         end
#     end
#
# end
# @show time() - t1
#
#
#
#
# return problem, x, x_k, P, q, err
# end
#
# problem, x, x_k, P, q, err  = solveit()
#
#
# @show x.value'*P*x.value + q'*x.value
#
# @show x_k'*P*x_k + q'*x_k


# mat"plot($err)"
