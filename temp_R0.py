import numpy as np
import matplotlib.pyplot as plt

# Constants
t = np.linspace(0, 365, 366)
T_0 = 5.73
om = 2*np.pi/365
epsilon = 1.91
fi = 92.3

fit = np.linspace(0, 100, 5)
a = np.linspace(100, 200, 101)
b = np.linspace(0.001, 0.3, 101)
beta_d1 = 2.14e-9
beta_i1 = 3.55e-9
Lambda = 2
d = 0.1/365
gamma1 = 0.14
gamma2 = 0.058

p1 = 10 ** 6.4
p2 = 10 ** 4

result1 = []
result2 = []
Temp = []

for i in t:
    # Temperature function
    temp = T_0 * (1 + epsilon * np.sin(om * i + fi))
    Temp.append(temp)
    
    Omega1 = 105.84 * np.exp(0.08 * temp)
    strand1 = ((beta_d1 * Lambda) / (d * (gamma1 + d))) + ((beta_i1 * Lambda * p1) / (d * Omega1 * (gamma1 + d)))
    
    Omega2 = 80 * np.exp(0.2 * temp)
    strand2 = ((beta_d1 * 100 * Lambda) / (d * (gamma2 + d))) + ((beta_i1 * 100 * Lambda * p2) / (d * Omega2 * (gamma2 + d)))
    
    result1.append(strand1)
    result2.append(strand2)

# Create a figure with 3 subplots (one for result1, one for result2, and one for temperature)
fig, ax = plt.subplots(1, 2, figsize=(10, 5), sharex=False)


# Plot the first result
ax[0].plot(Temp, result1, label='R0 of LPAI', color='b')
ax[0].set_ylabel('R0 value')
ax[0].set_xlabel('Temperature(°C)')


# Plot the second result

ax[0].plot(Temp, result2, label='R0 of HPAI', color='r')
ax[0].set_ylabel('R0 Value')
ax[0].set_xlim(np.min(Temp), np.max(Temp))
ax[0].grid(True)
ax[0].legend()
ax[0].axhline(y=1, color='orange', linestyle='--', label="threshold")  # Add horizontal line




# Plot the temperature
ax[1].plot(t, Temp, label='Temperature (°C)', color='g')
ax[1].set_xlabel('Days')
ax[1].set_ylabel('Temperature (°C)')
ax[1].set_xlim(0, 366)
ax[1].grid(True)
ax[1].legend()

# Adjust layout for better spacing
plt.tight_layout()
plt.savefig('R0_temp.png', dpi = 300)
plt.show()
