import numpy as np

# I-BVP
def f(x): # Initial value at t = 0
    return np.cos(np.pi * x / 2)

def a(t): # Boundary value at x = 0
    return 0

def b(x): # Boundary value at x = 1
    return 0

def ftcs(u, r, n, k): # FTCS scheme
    u[n + 1][k] = r * (u[n][k - 1] + u[n][k + 1]) + (1 - 2 * r) * u[n][k]

nu = 1.0 # Thermal diffusivity
x = 1
def sol(x, t): # Analytic solution
    return np.exp(-nu * np.pi ** 2 * t / 4) * np.cos(np.pi * x / 2)

n_pairs = int(input("Number of pairs (M, dt): "))
n_t = int(input("Number of time values t: "))
M_arr = [int(input("M = ")) for _ in range(n_pairs)]
dt_arr = [float(input("dt = ")) for _ in range(n_pairs)]
t_arr = [float(input("t = ")) for _ in range(n_t)]

for i in range(n_pairs):
    for j in range(n_t):
        M = M_arr[i]
        dt = dt_arr[i]
        t = t_arr[j]
        dx = float(x / M) # Space step size
        step_t = int(round(t / dt)) # Number of time steps
        r = nu * dt / (dx ** 2)
        u_first = [[0 for k in range(M + 1)] for n in range(step_t + 1)] # u[n][k]=u^n_k=u(k delta_x, n delta_t) approxes v(k delta_x,n delta_t)
        u_sec = [[0 for k in range(M + 1)] for n in range(step_t + 1)] # Second order Neumann approx

        for k in range(M + 1): # Initial value initialization
            u_first[0][k] = f(k * dx)
            u_sec[0][k] = f(k * dx)

        # Main loop (FTCS scheme)
        for n in range(step_t): # n = 0 1 ... step_t
            u_first[n + 1][M] = u_sec[n + 1][M] = b((n + 1) * dt) # Dirichlet right boundary
            for k in range(1, M): # k = 1 2 ... M - 1
                ftcs(u_first, r, n, k)
                ftcs(u_sec, r, n, k)
            u_first[n + 1][0] = u_first[n + 1][1] # First-order Neumann left boundary
            u_sec[n + 1][0] = u_sec[n][0] + 2 * r * (u_sec[n][1] - u_sec[n][0]) # Second-order Neumann left boundary

        # Error analysis
        v = np.array([sol(k * dx, t) for k in range(M + 1)]) # Analytic solution v^k(t), k = 0 1 ... M
        err_first = np.linalg.norm(v - np.array(u_first[step_t]), ord = np.inf) # Sup error first-order
        err_sec = np.linalg.norm(v - np.array(u_sec[step_t]), ord = np.inf) # Sup error second-order
        err_fs = np.linalg.norm(np.array(u_sec[step_t]) - np.array(u_first[step_t]), ord = np.inf)

        print("# (M, dt, t) = ", M, dt, t)
        print("# First-order error: ", err_first)
        print("# Second-order error: ", err_sec)
        print("# First-second error:", err_fs)
        print()

# (M, dt, t) =  20 0.001 0.06
# First-order error:  0.0168097592924501
# Second-order error:  9.198629266615743e-05
# First-second error: 0.016717772999783942

# (M, dt, t) =  20 0.001 0.1
# First-order error:  0.020003678065946495
# Second-order error:  0.00013889713139281223
# First-second error: 0.019864780934553683

# (M, dt, t) =  20 0.001 0.9
# First-order error:  0.01330771150878611
# Second-order error:  0.00017352579600697637
# First-second error: 0.013134185712779134

# (M, dt, t) =  40 0.00025 0.06
# First-order error:  0.00807735978635149
# Second-order error:  2.2976709995292666e-05
# First-second error: 0.008054383076356197

# (M, dt, t) =  40 0.00025 0.1
# First-order error:  0.009689772238488215
# Second-order error:  3.469521539889442e-05
# First-second error: 0.00965507702308932

# (M, dt, t) =  40 0.00025 0.9
# First-order error:  0.006683191988149109
# Second-order error:  4.336826160920848e-05
# First-second error: 0.006639823726539901