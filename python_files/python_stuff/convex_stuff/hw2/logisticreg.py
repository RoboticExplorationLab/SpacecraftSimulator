import cvxpy as cp
import numpy as np
import matplotlib.pyplot as plt
import math
import pdb

# load in problem data
DATA = np.loadtxt(
    '/Users/kevintracy/Documents/GitHub/python_stuff/convex_stuff/hw2/samples.csv',
    delimiter=',')

# get x's and y's
X = DATA[:, 0:-1]
Y = DATA[:, -1]

# create new matrix for convienence, b_i = -y_i*x_i
B = -np.diag(Y) @ X


# function value
def f(w, B):
    out = 0.0

    for i in range(np.size(B, 0)):
        b = B[i, :]
        t0 = math.exp(b.T @ w)

        out += math.log(1 + t0)

    return out


# function gradient
def df(w, B):
    out = np.zeros(10)

    for i in range(np.size(B, 0)):
        b = B[i, :]
        t0 = math.exp(b.T @ w)

        out += ((t0) / (1 + t0)) * b

    return out


# problem sizes
d = 10
N = 50

# solve with cvxpy
# w = cp.Variable(d)

# # loss function
# J = 0.0
# for i in range(N):
#     J += cp.logistic(B[i].T @ w)
# objective = cp.Minimize(J)

# # solve problem
# prob = cp.Problem(objective)

# result = prob.solve()

# print(result)

# print(w.value)

# print(f(w.value, B))
# print(df(w.value, B))

# now let's try GD with momentum


def gd_momentum(gamma, eta, B):
    w = np.zeros(d)

    g = np.zeros(d)
    f_vec = np.zeros(1001)

    # eta = 100
    # gamma = .01

    for i in range(1000):

        f_vec[i] = f(w, B)

        g = (1 - gamma) * g + gamma * df(w, B)

        w = w - eta * g

    f_vec[-1] = f(w, B)
    return f_vec / 50


gamma = .01

eta = .01
f_vec_1 = gd_momentum(gamma, eta, B)

eta = .1
f_vec_2 = gd_momentum(gamma, eta, B)

eta = 1.0
f_vec_3 = gd_momentum(gamma, eta, B)

eta = 10.0
f_vec_4 = gd_momentum(gamma, eta, B)

eta = 100.0
f_vec_5 = gd_momentum(gamma, eta, B)

# plt.semilogy(x, f_vec_1, label='\eta')
plt.figure(1)
plt.semilogy(f_vec_1, label='eta = 0.01')
plt.semilogy(f_vec_2, label='eta = 0.1')
plt.semilogy(f_vec_3, label='eta = 1')
plt.semilogy(f_vec_4, label='eta = 10')
plt.semilogy(f_vec_5, label='eta = 100')
plt.xlabel('t (iterations)')
plt.ylabel('f(w_t)')
plt.title("GD with Momentum (gamma = 0.01)")
plt.legend()

gamma = .1

eta = .01
f_vec_1 = gd_momentum(gamma, eta, B)

eta = .1
f_vec_2 = gd_momentum(gamma, eta, B)

eta = 1.0
f_vec_3 = gd_momentum(gamma, eta, B)

eta = 10.0
f_vec_4 = gd_momentum(gamma, eta, B)

# plt.semilogy(x, f_vec_1, label='\eta')
plt.figure(2)
plt.semilogy(f_vec_1, label='eta = 0.01')
plt.semilogy(f_vec_2, label='eta = 0.1')
plt.semilogy(f_vec_3, label='eta = 1')
plt.semilogy(f_vec_4, label='eta = 10')
plt.xlabel('t (iterations)')
plt.ylabel('f(w_t)')
plt.title("GD with Momentum (gamma = 0.1)")
plt.legend()

gamma = 1.0

eta = .01
f_vec_1 = gd_momentum(gamma, eta, B)

eta = .1
f_vec_2 = gd_momentum(gamma, eta, B)

eta = 1.0
f_vec_3 = gd_momentum(gamma, eta, B)

# plt.semilogy(x, f_vec_1, label='\eta')
plt.figure(3)
plt.semilogy(f_vec_1, label='eta = 0.01')
plt.semilogy(f_vec_2, label='eta = 0.1')
plt.semilogy(f_vec_3, label='eta = 1')
plt.xlabel('t (iterations)')
plt.ylabel('f(w_t)')
plt.title("GD with Momentum (gamma = 1)")
plt.legend()

plt.show()
