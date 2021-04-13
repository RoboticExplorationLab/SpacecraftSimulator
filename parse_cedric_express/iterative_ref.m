clear


n = 50;
A = randn(n,n);A = A'*A;

b = randn(n);
rho = 1e-4;
x1 = (A + 