import numpy as np
from scipy import integrate

# Constants
a = 0.01

# Time array
t = np.linspace(0, 365, 366)
omega = 2 * np.pi / 365

# Temperature function
def T(t):
    return 5.73 * (1 + 1.91 * np.sin(omega * t + 92.3))

# Parameters
gamma_1 = 0.14 
p_1 = 10**4
d1 = (0.1 / 365) + 0.5
d = (0.1 / 365)
Lambda = 2

# Transmission parameters
beta_d1 = 2.24e-9
beta_i1 = 3.55e-9


# Omega function
def Omega_1(T):
    return np.log(10) * np.exp(a* T - 3.5)

#np.log(10) * np.exp(abs(a) * T - abs(b))
# V matrix
def V_t(t):
    temp = T(t)
    omega1 = Omega_1(temp)
    return np.array([
        [gamma_1 + d1, 0],
        [-p_1, omega1]
    ], dtype=float)

# Constant F matrix
F_t = np.array([
    [beta_d1 * (Lambda / d),  beta_i1 * (Lambda / d)],
    [0, 0]
], dtype=float)

# Solver for ODE system
def dwdt(t, w, theta):
    V = V_t(t)
    dwdt = np.dot((-V + F_t / theta), w)
    return dwdt

# Solve with given initial condition and theta
def solver(w0, theta):
    sol = integrate.solve_ivp(
        lambda t, w: dwdt(t, w, theta),
        (0, t[-1] * 2),
        w0,
        method='BDF',
        max_step=1.0
    )
    return sol.y[:, -1]

# Identity matrix for initial condition
Y = np.eye(2)

# Sweep theta values
theta_values = np.linspace(2.5, 3, 500)  # Adjust range as needed
best_theta = None

for theta_current in theta_values:
    print(theta_current)
    Wmatrix = np.zeros([2, 2])
    for i in range(2):
        Wmatrix[:, i] = solver(Y[:, i], theta_current)

    # Compute eigenvalues and Emax
    eigenvalues = np.linalg.eigvals(Wmatrix)
    Emax = np.max(np.abs(eigenvalues))
    print('---', Emax)

    if 0.975 <= Emax <= 1.025:
        best_theta = theta_current
        print(f"Found theta: {best_theta} with Emax: {Emax}")
        break

