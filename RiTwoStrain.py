import numpy as np
from scipy.integrate import solve_ivp

t = np.linspace(0, 365, 366)
omega = 2*np.pi/365

def T(t):
    return 5.73 * (1 + 1.91 * np.sin(omega * t + 92.3))

gamma = 0.1
p_2 = 1e4
p_1 = 1e3
d1 = (0.1 / 365) + 0.88
d = (0.1 / 365)
Lambda = 2
alpha = 0.065

beta_d1 = 2.14e-9
beta_i1 = 3.55e-9

def Omega_2(temp):
    return np.log(10) * np.exp(0.0587 * temp - 3.6348)

def Omega_1(temp):
    return np.log(10) * np.exp(0.114*temp -  3.7594)

def V_t(t):
    temp = T(t)
    Omega1 = Omega_1(temp)
    Omega2 = Omega_2(temp)
    return np.array([
        [gamma + d, 0, 0, 0],
        [0, gamma + d1, 0, 0],
        [-p_1, 0, Omega1, 0],
        [0, -p_2, 0, Omega2]
    ], dtype=float)

Sstar = Lambda / d

F_t = np.array([
    [alpha* beta_d1 * Sstar,  (1-alpha)* beta_d1 * Sstar, alpha* beta_i1 * Sstar, (1-alpha)* beta_i1 * Sstar],
    [(1-alpha)* beta_d1 * Sstar, alpha* beta_d1 * Sstar, (1-alpha)* beta_i1 * Sstar, alpha* beta_i1 * Sstar],
    [0, 0, 0, 0],
    [0, 0, 0, 0]
], dtype=float)

def dwdt(t, w, theta):
    V = V_t(t)
    A = (-V + F_t / theta)
    return A @ w

def solver(w0, theta):
    sol = solve_ivp(
        lambda tt, ww: dwdt(tt, ww, theta),
        (0.0, 365.0),
        w0,
        method="BDF",
        max_step=1.0
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return sol.y[:, -1]

Y = np.eye(4)

def f(theta):
    W = np.zeros((4, 4))
    for i in range(4):
        w0 = Y[:, i]                 # use columns (basis vectors)
        W[:, i] = solver(w0, theta)  # column i of monodromy matrix

    eigvals = np.linalg.eigvals(W)
    spec_rad = np.max(np.abs(eigvals))   # spectral radius
    return spec_rad - 1.0

def find_theta(theta_low = 1, theta_high=5.0, tol=1e-6, max_expand=60, max_iter=80):
    f_low = f(theta_low)
    f_high = f(theta_high)

    # Expand upper bound until sign change
    k = 0
    while f_low * f_high > 0 and k < max_expand:
        theta_high *= 2.0
        f_high = f(theta_high)
        k += 1
    
    # Warning if there is no root 
    if f_low * f_high > 0:
        raise RuntimeError("Could not bracket root: f(theta_low) and f(theta_high) have same sign.")

    # Bisection method 
    for _ in range(max_iter):
        theta_mid = 0.5 * (theta_low + theta_high)
        f_mid = f(theta_mid)

        if abs(f_mid) < 1e-10 or (theta_high - theta_low) < tol:
            return theta_mid

        if f_mid > 0:
            theta_low = theta_mid
        else:
            theta_high = theta_mid

    return f_mid

theta_star = find_theta()
print("theta =", theta_star)
 
 