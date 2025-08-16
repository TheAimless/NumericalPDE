import numpy as np

def f(x): # Initial value at t = 0
    return np.sin(np.pi * x * 4)

def a(t): # Boundary value at x = 0
    return 0

def b(x): # Boundary value at x = 1
    return 0

nu = float(input("nu = ")) # Thermal diffusivity
def sol(x, t): # Analytic solution
    return np.exp(-16 * (np.pi ** 2) * nu * t) * np.sin(np.pi * x * 4)

delta_x = float(input("Δx = "))
delta_t = float(input("Δt = "))
x = 1
t = float(input("t = "))
M = int(1/delta_x) # Number of spatial steps
step_t = int(round(t / delta_t)) # Number of time steps
r = nu * delta_t / (delta_x ** 2)

u = [[0 for k in range(M+1)] for n in range(step_t + 1)] #u[n][k]=u^n_k=u(k delta_x, n delta_t) approxes v(k delta_x,n delta_t)

for k in range(M + 1): # Initial value initialization
    u[0][k] = f(k * delta_x)

# Main loop (FTCS scheme)
for n in range(step_t): # n = 0 1 ... step_t
    u[n + 1][0] = a((n + 1) * delta_t) # Dirichlet left boundary
    u[n + 1][M] = b((n + 1) * delta_t) # Dirichlet right boundary
    for k in range(1, M): # k = 1 2 ... M - 1
        u[n + 1][k] = (1 - 2 * r) * u[n][k] + r * (u[n][k + 1] + u[n][k - 1])

v = np.array([sol(k * delta_x, t) for k in range(M + 1)]) # Analytic solution v^k(t), k = 0 1 ... M
err = np.linalg.norm(v - np.array(u[step_t]), ord = np.inf) # Sup norm error

print("Error: ", err)

# Results:
# nu = .1:
# (i). dx = .1, dt = .05
  # t = .05: err = 1.379e-1
  # t = .1: err = 1.052e-1

# (ii). dx = .05, dt = .0125
# t = .05: err = 2.440e-2
# t = .1: err = 2.153e-2

# (i). dx = .01, dt = .0005
# t = .05: err = 9.447e-4
# t = .1: err = 8.569e-4