import numpy as np 


def Temp(t, T0, epsilon, om, fi):
    '''
    T = T0 + epsilon * np.sin(om * t + fi)
    ----------
    t : array or scalar
        time.
    T0 : scalar
        Average temperature.
    om : scalar
        Frequency.
    fi : scalar
        phase shift.
    epsilon : scalar
        Amplitude.

    Returns
    -------
    T : Instantaneous temperature value of temperature function T(t)
    '''
    T = T0 + epsilon*np.sin(om*t + fi)
    
    return T
    
       
    
def Viral(a, b, Temp):
    '''
    Input parameters
    ----------
    a : scalar
        Slope of viral decay.
    b : scalar
        Intercept of viral decay.
    Temp : Temp(t), scalar or array
        Obtained temperature via fitted temperature models.

    Returns
    -------
    Viral: Scalar or array
        Returned viral decay time at Temp. Depends on the input type 
        for parameter Temp, the return value could be the instantaneous decay rate 
        or a function over time. 

    '''  
    Viral = np.log(10)*np.exp(a*Temp - b)
    
    return Viral





def Avian(X, t, a, b, T0, epsilon, om, fi, beta_d = 2.13e-9, beta_i = 3.55e-9,
          d = 0.1/365, d2 = 0.88, p1 = 1e3, p2 = 1e4,
          gamma = 0.14, lambd = 2, eta = 0.038, alpha1 = 0.065, alpha2 = 0.065):
    '''
    Parameters
    ----------
    X : array
        Initial conditions of the dynamics.
    t : array
        Time period of the dynamics.
    beta_d : scalar
        Direct transmission rate. The default is 2.13e-9.
    beta_i : scalar
        Indirect transmission rate. The default is 3.55e-9.
    d : scalar
        Natura; death rate. The default is 0.1/365.
    d2 : scalar
        High pathogenic disease death rate. The default is 0.88.
    p1 : scalar
        Shedding rate for LPAI. The default is 1e3.
    p2 : scalar
        Shedding rate for HPAI. The default is 1e4.
    gamma : scalar
        Recovery rate. The default is 0.14.
    lambd : scalar
        Natural birth rate. The default is 2.
    eta : scalar
        Immunity waning rate. The default is 0.038.
    alpha1 : scalar
        Mutation rate for LPAI. The default is 0.065.
    alpha2 : scalar
        Mutation rate for HPAI. The default is 0.065.

    Returns
    -------
    Array dynamics
    '''
    current_temp = Temp(t, T0, epsilon, om, fi)
    S, I1, I2, R1, R2, I12, I21, R12, V1, V2 = X
    

    dotS = (lambd + eta*(R1 + R2 + R12)
            - beta_d*S*I1 - beta_d*I2*S
            - beta_i*V1*S - beta_i*V2*S
            - beta_d*S*(I12 + I21) - d*S)
    
    dotI1 = ((1- alpha1)*(beta_d*I1*S + beta_i*V1*S + beta_d*I21*S)
             + alpha2*(beta_d*I2*S + beta_d*V2*S + beta_d*S*I12)
             - gamma*I1 - d*I1)
    
    dotI2 = ((1- alpha2)*(beta_d*I2*S + beta_i*V2*S + beta_d*I12*S)
             + alpha1*(beta_d*I1*S + beta_d*V1*S + beta_d*S*I21)
             - gamma*I2 - (d+d2)*I2)
    
    dotR1 = (gamma*I1 - eta*R1
             - (1-alpha2)*R1*(V2*beta_i+ beta_d*I2 + beta_d*I12)
             - alpha1*R1*(beta_i*V1 + beta_d*I1 + beta_d*I21)
             - (d + d2)*R1)
    
    dotR2 = (gamma*I2 - eta*R2
             - (1-alpha1)*R2*(beta_i*V1 + beta_d*I1 + beta_d*I21)
             - alpha2*R2*(beta_i*V2 + beta_d*I2 +beta_d*I12)
             - d*R2)
    
    dotI12 = ((1-alpha2)*R1*(V2*beta_i+ beta_d*I2 + beta_d*I12)
              + alpha1*R1*(beta_i*V1 + beta_d*I1 + beta_d*I21)
              - (d + d2)*I12 - gamma*I12)
    
    dotI21 = ((1-alpha1)*R2*(beta_i*V1 + beta_d*I1 + beta_d*I21)
              + alpha2*R2*(beta_i*V2 + beta_d*I2 +beta_d*I12)
              - d*I21 - gamma*I21)
    
    dotR12 =  gamma*I12 + gamma*I21 - eta*R12 - d*R12
    
    dotV1 = p1*I1 + p1*I21 - Viral(0.114, 3.7594, current_temp) *V1  # regular LPAI
    dotV2 = p2*I2 + p2*I12 - Viral(a, b, current_temp)*V2  # HPAI
    
    return np.array([dotS, dotI1, dotI2, dotR1, dotR2, dotI12, dotI21, dotR12, dotV1, dotV2])


def DynamicsSolver(X0, t, a, b, T0, epsilon, om, fi, **kwargs):
    '''
    This function solves the ODE system and plots the dynamics
    Parameters
    ----------
    X0 : array
        All initial conditions in this exact order [S, I1, I2, R1, R2, I12, I21, R12, V1, V2].
    t : array
        Time period of the dynamics.
    a : scalar
        Slope of the viral decay rate.
    b : scalar
        Intercept of the viral decay rate.
    T0 : scalar
        Average temperature.
    epsilon : scalar
        Temperature amplitude.
    om : scalar
        Temperature frequency.
    fi : scalar
        Phase shift of temperature.
    **kwargs : optional
        Additional parameters to pass to avain().

    Returns
    -------
    Solution: array
        Solution array (shape len(t) x 10 ) and the dynamcis plot
    '''
    from scipy.integrate import odeint
    
    res = odeint(Avian, X0, t, args = (a, b, T0, epsilon, om, fi))
    S, I1, I2, R1, R2, I12, I21, R12, V1, V2 = res.T
    
    import matplotlib.pyplot as plt 
    
    #convert time into years
    years = t[-1]/365
    
    fig, ax = plt.subplots(figsize=(10, 8))
    plt.title(f"Limit cycle t = {years:.2f} years", fontsize = 20)
    plt.plot(t, S, label='Susceptible', color = 'black', linestyle = '-')
    plt.plot(t, I1,  color = 'b', linestyle = '--')
    plt.plot(t, I2,  color = 'r', linestyle = '--')
    plt.plot(t, I12,  color = 'tomato', linestyle = '-.')
    plt.plot(t, I21,  color = 'steelblue', linestyle = '-.')
    plt.plot(t, R1,  color = 'greenyellow', linestyle = ':')
    plt.plot(t, R2,  color = 'yellowgreen', linestyle = ':')
    plt.plot(t, R12, color = 'g', linestyle = ":")
    plt.xlabel('time (days)', fontsize = 20)
    plt.ylabel('Number of birds', fontsize = 20)
    plt.legend()
    plt.tight_layout()
    plt.show()

    fig, ax = plt.subplots(figsize=(10, 8))

    ax.set_title("Limit cycle - Infected and Recovered birds", fontsize=20)
    ax.plot(t, I1,  label='Infected by average LPAI', color='b',         linestyle='--')
    ax.plot(t, I2,  label='Infected by H5N1',         color='r',         linestyle='--')
    ax.plot(t, I12, label='Cross Infected by average LPAI', color='tomato',    linestyle='-.')
    ax.plot(t, I21, label='Cross Infected by H5N1',        color='steelblue', linestyle='-.')
    ax.plot(t, R1,  label='Recovered (LPAI)',              color='greenyellow', linestyle=':')
    ax.plot(t, R2,  label='Recovered (HPAI)',              color='yellowgreen', linestyle=':')
    ax.plot(t, R12, label='Recovered (cross)',             color='g',           linestyle=':')

    plt.xlabel('time (days)', fontsize = 20)
    plt.ylabel('Number of birds', fontsize = 20)
    plt.legend()
    plt.tight_layout()
    plt.show()


    fig, ax = plt.subplots(figsize=(10, 8))
    ax.plot(t, V1, label='LPAI particles in water', color='b')
    ax.plot(t, V2, label='H5N1 particles in water', color='r')

    plt.legend()
    plt.tight_layout()
    plt.show()

