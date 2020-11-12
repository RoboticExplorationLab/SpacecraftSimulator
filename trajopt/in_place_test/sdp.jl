using LinearAlgebra



function sdp_vec(Q,n)

    q = zeros(Int((n^2 - n)/2 + n))

    iter = 0
    for i = 1:n
        for j = 1:i
            iter+=1
            q[iter] = Q[j,i]
        end
    end
    return q
end

function sdp_mat(q,n)

    Q = zeros(n,n)

    iter = 0
    for i = 1:n
        for j = 1:i
            iter+=1
            Q[j,i] = q[iter]
        end
    end
    return (Q + Q') - diagm(diag(Q))
end

function sdp_proj(q,n)
    Q = sdp_mat(q,n)

    V,S = eigen(Q)

    return S*diagm(max.(V,0))*S'

end


n = 3

Q = randn(n,n);Q = Q'*Q


q = sdp_vec(Q,n)

Q2 = sdp_mat(q,n)


Q = randn(5,5)
Q = .5*(Q + Q')

q = sdp_vec(Q,5)

Qn = sdp_proj(q,5)
