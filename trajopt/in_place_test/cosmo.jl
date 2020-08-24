using LinearAlgebra, Convex, SCS, SparseArrays, Infiltrator, QDLDL


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

function cone_proj(x::Vector)::Vector
    """Project onto the SOC"""


    v = x[1:(end-1)]
    s = x[end]

    nv = norm(x[1:end-1])
    if nv<=-x[end]
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
    fA2 = qdldl(QCQP.LS_A)
    nu = length(q)

    for k = 1:30


        # solve linear system
        QCQP.LS_sol .= solve(fA2,[(-QCQP.q + QCQP.σ*QCQP.x_k);(QCQP.b-QCQP.s_k + (1/QCQP.ρ)*QCQP.y_k)])

        # put the solution to the linear system in the correct spots
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

        # termination criteria
        if rem(k,5)==0
            if norm(QCQP.A*QCQP.x_k + QCQP.s_k - QCQP.b)<1e-4
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

    t1 = time()
    Convex.solve!(problem, SCS.Optimizer)
    convex_time1 = time()-t1
    # @infiltrate
    # error()

    @show x.value

    # my part of it
    t1 = time()
    QCQP_solve!(QCQP::QCQP_struct,P,q)
    my_time1 = time()-t1

    @show QCQP.x_k


    #now that it is warm, let's speed test it
    N_trials = 1000
    convex_times = zeros(N_trials)
    my_times = zeros(N_trials)

    # for i = 1:N_trials
    #
    #     P_new = .1*randn(3,3)
    #     P += P_new'*P_new
    #
    #     # solve again
    #     t1 = time()
    #     Convex.solve!(problem, SCS.Optimizer, warmstart=true)
    #     convex_times[i] = time()-t1
    #     # @infiltrate
    #     # error()
    #
    #     @show x.value
    #
    #     # my part of it
    #     t1 = time()
    #     QCQP_solve!(QCQP::QCQP_struct,P,q)
    #     my_times[i] = time()-t1
    # end
    t1 = time()
    P_first = copy(P)
    for i = 1:N_trials

        P_new = .1*randn(3,3)
        P += P_new'*P_new

        # solve again
        # t1 = time()
        Convex.solve!(problem, () -> SCS.Optimizer(verbose=false), warmstart=true)
        # convex_times[i] = time()-t1
    end
    convex_times = time() - t1
    P = copy(P_first)
    t2 = time()
    for i = 1:N_trials
        P_new = .1*randn(3,3)
        P += P_new'*P_new
        # @infiltrate
        # error()

        # @show x.value

        # my part of it
        # t1 = time()
        QCQP_solve!(QCQP::QCQP_struct,P,q)
        # my_times[i] = time()-t1
    end
    my_times = time() - t2

    return convex_times, my_times
end


# mat"
# figure
# hold on
# plot($convex_times)
# plot($my_times)
# legend('Convex','Mine')
# hold off
# end"
