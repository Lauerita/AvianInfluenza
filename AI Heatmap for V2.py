import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

def avian(X, t):
    S, I1, I2, R1, R2, V1, V2 = X
    Temp = T_0 * (1 + epsilon * np.sin(om*t + fi))
    dotS = lambd*S + eta*(R1 + R2) - S*(beta_d_1*I1 + beta_i_1*V1 + beta_d_2*I2 + beta_i_2*V2 + d1)
    dotI1 = S*(beta_d_1*I1 + beta_i_1*V1) + R2*(beta_d_1*I1 + beta_i_1*V1) - I1*(gamma1)
    dotI2 = S*(beta_d_2*I2 + beta_i_2*V2) + R1*(beta_d_2*I2 + beta_i_2*V2) - I2*(gamma2 + d2)
    dotR1 = I1*gamma1 - R1*(beta_d_2*I2 + beta_i_2*V2 + d1 + eta)
    dotR2 = I2*gamma2 - R2*(beta_d_1*I1 + beta_i_1*V1 + d1 + eta)
    dotV1 = p1*I1 - a1*np.exp(-b1*Temp)*V1 
    dotV2 = p2*I2 - a2*np.exp(-b2*Temp)*V2
    return np.array([dotS, dotI1, dotI2, dotR1, dotR2, dotV1, dotV2])

# Parameter setup
grid = 20
a = np.linspace(50, 200, grid)
b = np.linspace(0.01, 0.1, grid)
storage = np.zeros([grid, grid])

# Model parameters
days = 365
t = np.linspace(1, days, days + 1)
N = 5000
I_01, I_02 = 200, 20
S_0 = N - (I_01 + I_02)
R_01 = R_02 = 0
V_01, V_02 = 10**4.7, 0

# Fixed parameters
lambd = 0.1/365
eta = 0.038
beta_d_1 = 2.14e-9
beta_i_1 = 3.55e-9
beta_d_2 = beta_d_1 * 46.55
beta_i_2 = beta_i_1 * 46.55
d1 = d2 = 0.1/365
gamma1 = 0.14
gamma2 = 0.058
p1 = 1*10**6.4 
p2 = 1e4
T_0 = 5.73
om = 2*np.pi/365
epsilon = 1.91
fi = 92.3
a1 = 105.84
b1 = 0.08

# Calculate values
for i in range(len(a)):
    for j in range(len(b)):
        a2, b2 = a[i], b[j]
        X_0 = [S_0, I_01, I_02, R_01, R_02, V_01, V_02]
        res = integrate.odeint(avian, X_0, t)
        S, I1, I2, R1, R2, V1, V2 = res.T
        # Convert to percentage (0-100 scale)
        storage[i, j] = 100 * I2[-1]/(I1[-1] + I2[-1])

# Create the plot
fig, ax = plt.subplots(figsize=(10, 8))

# Plot heatmap with updated scale
im = plt.imshow(storage, interpolation='spline36', cmap='coolwarm', 
                vmin=0, vmax=100)  # Set explicit range for colormap

# Add contour lines (adjusted for percentage scale)
contour = ax.contour(storage, levels=[40, 50, 60], colors='black', linewidths=1)
ax.clabel(contour, inline=True, fontsize=8)

# Customize axes
plt.xticks(np.arange(0, grid, 2), labels=[f'{round(i, 1)}' for i in a[::2]], 
           rotation=45, fontsize=8)
plt.yticks(np.arange(0, grid, 2), labels=[f'{round(i, 3)}' for i in b[::2]], 
           fontsize=8)

# Add colorbar with percentage labels
cbar = plt.colorbar(im)
cbar.set_label('Percentage of HPAI (%)', fontsize=10)

# Labels and title
plt.xlabel('a-values', fontsize=10)
plt.ylabel('b-values', fontsize=10)
plt.title('Percentage of HPAI in the water', fontsize=12)

# Adjust layout and save
plt.tight_layout()
plt.savefig('heatmap_with_line.png', dpi=300, bbox_inches='tight')

