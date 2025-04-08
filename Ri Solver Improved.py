import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

# Time array
t = np.linspace(0, 365, 366)
omega = 2 * np.pi / 365

# Temperature function
def T(t, a = 5.73):
    return a * (1 + 1.91 * np.sin(omega * t + 92.3))

# Define parameters
gamma_1 = 0.14
gamma_2 = 0.058
p_1 = 10 ** 6.4
p_2 = 1e4
d = 0.1 / 365
Lambda = 2
beta_d1 = 2.24e-9
beta_i1 = 3.55e-9
beta_d2 = 2.24e-7
beta_i2 = 3.55e-7

# make theta into an array 
theta = np.linspace(1, 10, 1000)

# Define time-dependent Omega functions
def Omega_1(T):
    return 105.84 * np.exp(0.08 * T)

def Omega_2(T):
    return 80 * np.exp(0.2 * T)

# Define the V matrix function
def V_t(t):
    temp = T(t)
    omega1 = Omega_1(temp)
    omega2 = Omega_2(temp)
    
    return np.array([
        [gamma_1 + d, 0, 0, 0, 0, 0],
        [0, gamma_2 + d, 0, 0, 0, 0],
        [0, 0, gamma_2 + d, 0, 0, 0],
        [0, 0, 0, gamma_1 + d, 0, 0],
        [-p_1, 0, 0, -p_1, omega1, 0],
        [0, -p_2, -p_2, 0, 0, omega2]
    ], dtype=float)

# Define F matrix (constant)
F_t = np.array([
    [beta_d1 * (Lambda / d), 0, 0, beta_d1 * (Lambda / d), beta_i1 * (Lambda / d), 0],
    [0, beta_d2 * (Lambda / d), beta_d2 * (Lambda / d), 0, 0, beta_i2 * (Lambda / d)],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 0, 0]
], dtype=float)

# Define the system of ODEs
def dwdt(t, w):
    V = V_t(t)
    F = F_t
    dwdt = np.dot((-V + F/theta_current), w)
    return dwdt

  
# define the solver for the differential equation
def solver(w0):
    sol = integrate.solve_ivp(
        dwdt, 
        (t[0], t[-1]), 
        w0, 
        method='BDF',  
        max_step=1.0
    )
    return sol.y[:, -1]


for theta_current in theta:
    #print('The current theta value is:', theta_current)
    
    
    # Initial condition identity matrix 
    Y = np.eye(6)
    
    #set up for loop for each column vector
    Wmatrix = np.zeros([6,6])

    for i in range(6):
        Wmatrix[:, i] = solver(Y[:, i])
        
        
    # Calculate all eigenvalues
    eigenvalues = np.linalg.eigvals(Wmatrix)

    # Find the maximum eigenvalue (by magnitude)
    Emax = np.max(np.abs(eigenvalues))
    #print('The current Emax is:', Emax)
    
    if 0.95 <= Emax <= 1.05:
        best_theta = theta_current
        print(f"Found theta: {best_theta} with Emax: {Emax}")
        break 
        
    




