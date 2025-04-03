import numpy as np
import matplotlib.pyplot as plt

# Constants
T_0 = 5.73
om = 2*np.pi/365
epsilon = 1.91
fi = 92.3
Temp = T_0 * (1 + epsilon * np.sin(om + fi))

fit = np.linspace(0, 100, 5)
a = np.linspace(50, 200, 101)
b = np.linspace(0.001, 0.3, 101)
beta_d1 = 2.14e-9
beta_i1 = 3.55e-9
Lambda = 2
d = 0.1 / 365
gamma1 = 0.14
gamma2 = 0.058
Omega1 = 105.84 * np.exp(0.08 * Temp)
p1 = 10 ** 3
p2 = 10 ** 4

# Create a list to store all the matrices for consistent color scali

# Create a single figure with subplots for each 'fit' value
fig, ax = plt.subplots(1, len(fit), figsize=(25, 8))

# Create a heatmap for each fit value
for idx, i in enumerate(fit):
    # Calculate beta values for current fit value
    beta_d2 = i * beta_d1
    beta_i2 = i * beta_i1
    
    # Initialize matrix for Strand 2
    matrix2 = np.zeros((len(a), len(b)))
    
    # Calculate matrix values
    for j in range(len(a)):
        for k in range(len(b)):
            # Calculate Omega2 for each a,b combination
            Omega2 = a[j] * np.exp(b[k] * Temp)
            
            # Calculate strand values for Strand 2
            strand1 = (beta_d1 * Lambda) / (d * (gamma1 + d)) + (beta_i1 * Lambda * p1) / (d * Omega1 * (gamma1 + d))
            strand2 = (beta_d2 * Lambda) / (d * (gamma2 + d)) + (beta_i2 * Lambda * p2) / (d * Omega2 * (gamma2 + d))
            
            # Store the strand values in their respective matrix
            matrix2[j, k] = strand2

    
    # Create the heatmap for Strand 2 with the same color scale across all subplots
    im2 = ax[idx].imshow(matrix2, origin='lower', aspect='auto', cmap = 'coolwarm',
                         extent=[b.min(), b.max(), a.min(), a.max()],
                         vmin=0, vmax=20)
    ax[idx].set_xlabel('b values')
    ax[idx].set_ylabel('a values')
    ax[idx].set_title(f'Strand 2 (fit_value = {i:.2f})')
    
    # Find contour for R0 = 7.69
    contour_levels = [6.25]
    contour = ax[idx].contour(b, a, matrix2, levels=contour_levels, colors='black', linewidths=2)
    ax[idx].clabel(contour, inline=True, fontsize=10, fmt='%.2f')


# Create a single colorbar for all subplots, with consistent scale
cbar = fig.colorbar(im2, ax=ax, orientation='horizontal', pad=0.5, shrink=0.5)
cbar.set_label('HPAI Reproduction Number')

# Adjust layout and show the plot
plt.tight_layout()
plt.suptitle(f'Reproduction Number of HPAI', fontsize=16, y=0.98)
plt.subplots_adjust(top=0.9, bottom=0.3)  # Adjust 'bottom' to create space for the colorbar
plt.savefig('R0_HPAI.png', dpi = 300)
plt.show()

