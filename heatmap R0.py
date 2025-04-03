import numpy as np
import matplotlib.pyplot as plt

# Constants
average = 20
fit = np.linspace(0, 100, 5)
a = np.linspace(150, 200, 101)
b = np.linspace(0.001, 0.3, 101)
beta_d1 = 2.14e-9
beta_i1 = 3.55e-9
Lambda = 2
d = 0.1/365
gamma1 = 0.14
gamma2 = 0.058
Omega1 = 105.84 * np.exp(-0.08 * average)
p1 = 10 ** 6.4
p2 = 10 ** 4

fig, ax = plt.subplots(1, figsize=(8, 8))
   
# Create a separate figure for each fit value
for i in fit:
    fig, ax = plt.subplots(1, figsize=(8, 8))
    # Create a new figure for each fit value 
    # Calculate beta values for current fit value
    beta_d2 = i * beta_d1
    beta_i2 = i * beta_i1
    
    # Initialize matrix
    matrix = np.zeros((len(a), len(b)))
    matrix2 = np.zeros((len(a), len(b)))
    
    # Calculate matrix values
    for j in range(len(a)):
        for k in range(len(b)):
            # Calculate Omega2 for each a,b combination
            Omega2 = a[j] * np.exp(-b[k] * average)
            
            # Calculate strand values
            strand1 = (beta_d1 * Lambda) / (d * (gamma1 + d)) + (beta_i1 * Lambda * p1) / (d * Omega1 * (gamma1 + d))
            strand2 = (beta_d2 * Lambda) / (d * (gamma2 + d)) + (beta_i2 * Lambda * p2) / (d * Omega2 * (gamma2 + d))
            
            # Store the strand values in their respective matrices
            matrix[j, k] = strand1
            matrix2[j, k] = strand2
    
    # Create the first heatmap
    #im1 = ax[0].imshow(matrix, origin='lower', aspect='auto', 
     #             extent=[b.min(), b.max(), a.min(), a.max()])
    #fig.colorbar(im1, ax=ax[0], label='Strand 1 Value')
    #ax[0].set_xlabel('b values')
    #ax[0].set_ylabel('a values')
    #ax[0].set_title(f'Strand 1 (fit_value = {i:.2f})')
    
    # Create the second heatmap
    im2 = ax.imshow(matrix2, origin='lower', aspect='auto', 
                  extent=[b.min(), b.max(), a.min(), a.max()])
    fig.colorbar(im2, ax=ax, label='Strand 2 Value')
    ax.set_xlabel('b values')
    ax.set_ylabel('a values')
    ax.set_title(f'Strand 2 (fit_value = {i:.2f})')
    
    plt.tight_layout()
    plt.suptitle(f'Strand Values Comparison (fit_value = {i:.2f})', fontsize=16, y=0.98)
    plt.subplots_adjust(top=0.9)
    plt.show()  # Display the current figure before moving to the next one