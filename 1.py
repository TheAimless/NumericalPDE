# HW 1.4.1
# Second-order scheme
import numpy as np

def f(x): # Initial value at t = 0
    return np.cos(np.pi * x / 2)

def a(x): # Boundary value at x = 0
    return 0

def b(x): # Boundary value at x = 1
    return 0

def sol(x, t): # Analytic solution
    return np.exp(-nu * (np.pi / 2) ** 2 * t) * np.cos(np.pi * x / 2)

def ftcs(u, r, n, k): # FTCS scheme
    u[n + 1][k] = r * (u[n][k - 1] + u[n][k + 1]) + (1 - 2 * r) * u[n][k]

nu = 1.0 # Thermal diffusivity
M = 20 # Number of spatial steps
t = float(input("t = "))
delta_t = float(input("delta_t = "))
x = 1
delta_x = 1/M
step_t = int(round(t / delta_t))
r = nu * (delta_t / delta_x ** 2)

u = [[0 for i in range(M+1)] for j in range(step_t + 1)] #u(n,k)=u^n_k approxes v(k delta_x,n delta_t)

for k in range(M + 1):
    u[0][k] = f(k * delta_x)

for n in range(step_t):
    u[n + 1][0] = u[n][0] + 2 * r * (u[n][1] - u[n][0]) # O(delta_x^2)
    u[n + 1][M] = b((n + 1)* delta_t)
    ftcs(u, r, n, k)

err = 0

for n in range(step_t + 1):
    for k in range(M + 1):
        u[n][k] = np.fabs(u[n][k]-sol(k * delta_x, n * delta_t))
        err = max(err, u[n][k])

for i in u:
    print(i)

print(err)