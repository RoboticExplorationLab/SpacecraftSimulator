clear


P = randn(2,2);P = P'*P;
c = 3.0;

cvx_begin 
    variable x(2)
    minimize(norm(x))
    subject to 
      