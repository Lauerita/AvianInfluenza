import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Time array
t = np.linspace(0, 365, 366)
omega = 2 * np.pi / 365

# Temperature function
def T(t):
    return 5.73 * (1 + 1.91* np.sin(omega * t + 92.3))

# Define parameters
gamma_1 = 0.14
gamma_2 = 0.058
p_1 = 10 ** 6.4
p_2 = 1e4
d = 0.1 / 365

# make theta into an array 
theta = 4.626
Lambda = 2

# Define transmission parameters
beta_d1 = 2.24e-9
beta_i1 = 3.55e-9
beta_d2 = 2.24e-7
beta_i2 = 3.55e-7

# Define time-dependent Omega functions
def Omega_1(T):
    return 105.84 * np.exp(0.08 * T)

def Omega_2(T):
    return 80 * np.exp(0.2 * T)

# Define the V matrix function
def V_t(t):
    temp = T(t)
    #omega1 = Omega_1(temp)
    omega2 = Omega_2(temp)
    
    return np.array([
        [gamma_2 + d, 0, 0],
        [0, gamma_2 + d, 0],
        [-p_2, -p_2, omega2]
    ], dtype=float)

# Define F matrix (constant)
F_t = np.array([
    [beta_d2 * (Lambda / d), beta_d2 * (Lambda / d), beta_i2 * (Lambda / d)],
    [0, 0, 0],
    [0, 0, 0]
], dtype=float)

# Define the system of ODEs
def dwdt(t, w):
    V = V_t(t)
    F = F_t
    dwdt = np.dot((-V + F/theta), w)
    return dwdt

# Initial condition identity matrix 
Y = np.eye(3)
    
    
# define the solver for the differential equation
def solver(w0):
    sol = integrate.solve_ivp(
        dwdt, 
        (0, t[-1]*2), 
        w0, 
        method='BDF',  # Use BDF method for stiff problems,
        max_step=1.0
    )
    return sol.y[:, -1]

#set up for loop for each column vector
Wmatrix = np.zeros([3,3])

for i in range(3):
    Wmatrix[:, i] = solver(Y[:, i])
    
    
# Calculate all eigenvalues
eigenvalues = np.linalg.eigvals(Wmatrix)

# Find the maximum eigenvalue (by magnitude)
Emax = np.max(np.abs(eigenvalues))
print(Emax)
