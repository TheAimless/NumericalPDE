import math

def f(x):
    return math.cos(math.pi * x / 2)

def a(x):
    return 0

def b(x):
    return 0

def sol(x, t):
    return math.exp(-nu * (math.pi / 2) ** 2 * t) * math.cos(math.pi * x / 2)

M = 10; nu = 1.0
t = 0.06; x = 1
delta_x = 1/M; delta_t = .004
step_t = int(round(t / delta_t))

u = [[0 for i in range(M+1)] for j in range(step_t + 1)] #u(n,k)=u^n_k approxes v(k delta_x,n delta_t)

for k in range(M + 1):
    u[0][k] = f(k * delta_x)

for n in range(step_t):
    u[n + 1][0] = u[n][0] + 2 * nu * delta_t / (delta_x ** 2) * (u[n][1] - u[n][0]) # O(delta_x^2)
    u[n + 1][M] = b((n + 1)* delta_t)
    for k in range(1, M):
        u[n + 1][k] = u[n][k] + nu * (delta_t / delta_x ** 2) * (u[n][k + 1] - 2 * u[n][k] + u[n][k - 1])
    #u[n + 1][0] = u[n + 1][1] # O(delta_x)

err = 0

for n in range(step_t + 1):
    for k in range(M + 1):
        u[n][k] = math.fabs(u[n][k]-sol(k * delta_x, n * delta_t))
        err = max(err, u[n][k])

for i in u:
    print(i)

print(err)