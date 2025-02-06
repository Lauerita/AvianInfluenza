import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import numba
from scipy import integrate


#----------Initial Conditions-------------#
days = 365*2
t = np.linspace(1, days, days + 1) # every day for 6 months 

N = 5000
I_01 = 200
I_02 = 250
S_0 = N - (I_01 + I_02)
R_01 = 0
R_02 = 0
V_01 = 10**4.7 # minimal viral load of low pathogenic virus 
V_02 = 0


#--------- Parameters-----------#

# Relative fitness Parameter 
#param = [223.1, 445.222, 667.333, 1333.666, 2000]
zeta = 41

#-----------#

lambd = 0.1/365
eta = 0.038

beta_d_1 = 2.14*10**-9
beta_i_1 = 3.55*10**-9
beta_d_2 = beta_d_1 * zeta## relative fitness to low pathogenic virus 0 <c< infinity 
beta_i_2 = beta_i_1 * zeta ## relative fitness to low pathogenic 

d1 = 0.1/365 # natural death rate per day 
d2 = d1 # high pathogenic death rate 

gamma1 = 0.14 # recovery rate of low pathogenic virus 
gamma2 = 0.058 #recovery rate of high pathogenic virus

p1 = 1*10**6.4 # shedding rate of low pathogenic virus
p2 = 1*10**4 # shedding rate of high pathogenic virus 


# Temperature parameters 
T_0 = 5.73
om = 2*np.pi/365
epsilon = 1.91
fi = 92.3

a1 = 105.84
b1 = -0.08

a2 = 120
b2 = -0.03

#T_may = T_0* (1 + epsilon* np.sin(om*May + fi))

#-------------------------------#
#
#
#
#---------- ODE system ----------------#

def avian(X, t):
    S, I1, I2, R1, R2, V1, V2 = X
    Temp = T_0* (1 + epsilon* np.sin(om*t + fi))
    dotS = lambd*S + eta*(R1 + R2) - S*(beta_d_1*I1 + beta_i_1*V1 + beta_d_2*I2 + beta_i_2*V2 + d1)
    dotI1 = S*(beta_d_1*I1 + beta_i_1*V1) + R2*(beta_d_1*I1 + beta_i_1*V1) - I1*(gamma1)
    dotI2 = S*( beta_d_2*I2 + beta_i_2*V2) + R1*( beta_d_2*I2 + beta_i_2 *V2) - I2*(gamma2 + d2)
    dotR1 = I1*gamma1 - R1*(beta_d_2*I2 + beta_i_2*V2 + d1 + eta)
    dotR2 = I2*gamma2 - R2*(beta_d_1*I1 + beta_i_1*V1 + d1 + eta)
    dotV1 = p1*I1 - a1*np.exp(b2*Temp)*V1 
    dotV2 = p2*I2 -  a2*np.exp(b2*Temp)*V2 #need to be different

    return  np.array([dotS, dotI1, dotI2, dotR1, dotR2, dotV1, dotV2])
    
    
 #---------------ODE solver-----------------------#


X_0 = S_0, I_01, I_02, R_01, R_02, V_01, V_02

res = integrate.odeint(avian, X_0, t)

S, I1, I2, R1, R2, V1, V2 = res.T


#--------------------GRAPHS------------------#

Temp = T_0* (1 + epsilon* np.sin(om*t + fi))
plt.figure()
plt.grid()
plt.title('Temperature Change in Helsinki, Finland')
plt.plot(t, Temp, label = 'temperature in the breeding ground')
plt.xlabel('Time t, [days]')
plt.ylabel('Temperature in Celsius')
plt.savefig('TempHelsinki.png', dpi = 300)

#----------------------#

plt.figure(figsize=(10, 8))
plt.grid()
plt.title("General Graphical Solution")
plt.plot(t, S, label='Susceptible')
plt.plot(t, I1, label='Infected by LPAI')
plt.plot(t, I2, label='Infected by HPAI')
plt.plot(t, R1, label='Recovered with immunity of LPAI')
plt.plot(t, R2, label='Recovered with immunity of HPAI')
plt.xlabel('Time t, [days]')
plt.ylabel('Numbers of birds')
plt.legend()
plt.savefig('GensolutionAI.png', dpi = 300)

#-----------------------#

plt.figure(figsize=(10, 8))
plt.grid()
plt.title("Infected birds and Recovered birds")
plt.plot(t, I1, label='Infected by LPAI', color = 'orange')
plt.plot(t, I2, label='Infected by HPAI', color = 'green')
plt.plot(t, R1, label='Recovered with immunity of LPAI', color = 'red')
plt.plot(t, R2, label='Recovered with immunity of HPAI', color = 'purple')
plt.xlabel('Time t, [days]')
plt.ylabel('Numbers of birds')
plt.xlim([0, days])
#plt.ylim(0, 100)
plt.legend()

#-------------------------#

fig, ax = plt.subplots(2, 1, figsize=(10, 8))

ax[0].plot(t, V1, label='LPAI particles in water', color = 'b')
ax[0].set_xlabel('Time t, [days]')
ax[0].set_ylabel('viral particle load')
ax[0].legend()

ax[1].plot(t, V2, label='HPAI particles in water', color = 'r')
ax[1].set_xlabel('Time t, [days]')
ax[1].set_ylabel('viral particle load')
ax[1].legend()

plt.tight_layout()  # Adjusts spacing between subplots
plt.savefig('ViralPart.png', dpi=300)
plt.show()



   