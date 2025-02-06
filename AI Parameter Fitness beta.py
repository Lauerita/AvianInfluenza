import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numba
from scipy import integrate

def avian(X, t, params):
    """ODE solver with parameters passed as argument"""
    S, I1, I2, R1, R2, V1, V2 = X
    (lambd, eta, beta_d_1, beta_i_1, beta_d_2, beta_i_2, d1, d2, 
     gamma1, gamma2, p1, p2, T_0, om, epsilon, fi) = params
    
    Temp = T_0 * (1 + epsilon * np.sin(om*t + fi))
    dotS = lambd*S + eta*(R1 + R2) - S*(beta_d_1*I1 + beta_i_1*V1 + beta_d_2*I2 + beta_i_2*V2 + d1)
    dotI1 = S*(beta_d_1*I1 + beta_i_1*V1) + R2*(beta_d_1*I1 + beta_i_1*V1) - I1*(gamma1)
    dotI2 = S*(beta_d_2*I2 + beta_i_2*V2) + R1*(beta_d_2*I2 + beta_i_2*V2) - I2*(gamma2 + d2)
    dotR1 = I1*gamma1 - R1*(beta_d_2*I2 + beta_i_2*V2 + d1 + eta)
    dotR2 = I2*gamma2 - R2*(beta_d_1*I1 + beta_i_1*V1 + d1 + eta)
    dotV1 = p1*I1 - 105.84*np.exp(-0.08*Temp)*V1 
    dotV2 = p2*I2 - 105.84*np.exp(-0.08*Temp)*V2
    
    return np.array([dotS, dotI1, dotI2, dotR1, dotR2, dotV1, dotV2])

# Set up all constant parameters
days = 365
t = np.linspace(1, days, days + 1)

# Population parameters
N = 5000
I_01 = 200
I_02 = 20
S_0 = N - (I_01 + I_02)
R_01 = 0
R_02 = 0
V_01 = 10**(4.7)
V_02 = 0

# Disease parameters
lambd = 0.1/365
eta = 0.038
beta_d_1 = 2.14e-9
beta_i_1 = 3.55e-9
d1 = 0.1/365
d2 = d1
gamma1 = 0.14
gamma2 = 0.058
p1 = 1*10**6.4
p2 = 1e4

# Temperature parameters
T_0 = 5.73
om = 2*np.pi/365
epsilon = 1.91
fi = 92.3

# Initial conditions
X_0 = [S_0, I_01, I_02, R_01, R_02, V_01, V_02]

# Parameter sweep setup
param = np.linspace(0, 100, 1000)
I1_sum = []
I2_sum = []

# Main loop
for i in param:
    # Update fitness-dependent parameters
    beta_d_2 = beta_d_1 * i
    beta_i_2 = beta_i_1 * i
    
    # Pack parameters for ODE solver
    params = (lambd, eta, beta_d_1, beta_i_1, beta_d_2, beta_i_2, d1, d2,
             gamma1, gamma2, p1, p2, T_0, om, epsilon, fi)
    
    # Solve ODE system
    res = integrate.odeint(avian, X_0, t, args=(params,))
    S, I1, I2, R1, R2, V1, V2 = res.T
    
    # Store results
    I1_sum.append(I1[-1])
    I2_sum.append(I2[-1])
    
    # Print when populations are close
    if abs(I1[-1] - I2[-2]) < 1:
        print(f"Equilibrium found at fitness value: {i:.2f}")

# Plotting
fig, ax = plt.subplots(figsize=(10, 8))
plt.grid(True)
plt.plot(param, I1_sum, label='LPAI')
plt.plot(param, I2_sum, label='HPAI')
plt.xlabel('Fitness value [0,100]')
plt.ylabel('Number of infected birds at the end of the season')
plt.legend()

# Save the figure
plt.savefig('ParameterSweepBeta.png', dpi=300, bbox_inches='tight')