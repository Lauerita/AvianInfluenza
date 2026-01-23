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
    epsilon : sclalr
        Amplitude.

    Returns
    -------
    T : Instantaenous temperature value of temperature function T(t)
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
    Viral: Scalalr or array
        Returned viral decay time at Temp. Depends on the input type 
        for parameter Temp, return value could be instantaneous decay rate 
        or a function over time. 

    '''  
    Viral = np.log(10)*np.exp(a*Temp - b)
    
    return Viral

   
# Define A and B
def R0(t, a, b, T0, epsilon, fi, om, beta_d = 2.3e-9, beta_i = 3.55e-9, p1 = 1e3, 
       p2 = 1e4, gamma = 0.24, d = 0.1/365, d2 = 0.88,
       alpha1 = 0.065, alpha2 = 0.065):
    '''
    Parameters
    ----------
    t : array or scalar
    Time points for simulation (days).
    a : scalar
    Slope of viral decay for HPAI.
    b : scalar
    Intercept of viral decay for HPAI.
    T0 : scalar
    Average temperature (°C).
    epsilon : scalar
    Temperature amplitude (°C).
    fi : scalar
    Phase shift for temperature function (radians).
    om : scalar
    Frequency of temperature oscillation (radians/day).
    beta_d : scalar, optional
    Direct transmission rate. The default is 2.3e-9.
    beta_i : scalar, optional
    Indirect transmission rate. The default is 3.55e-9.
    p1 : scalar, optional
    Viral shedding rate for LPAI. The default is 1e3.
    p2 : scalar, optional
    Viral shedding rate for HPAI. The default is 1e4.
    gamma : scalar, optional
    Recovery rate. The default is 0.24.
    d : scalar, optional
    Natural death rate. The default is 0.1/365.
    d2 : scalar, optional
    Disease-induced death rate for HPAI. The default is 0.88.
    alpha1 : scalar, optional
    Mutation rate from LPAI to HPAI. The default is 0.065.
    alpha2 : scalar, optional
    Mutation rate from HPAI to LPAI. The default is 0.065.

    Returns
    -------
    R0 : scalar
        The time-invariant basic reproduction number for both viruses 

    '''
    current_temp = Temp(t, T0, epsilon, om, fi)
    S_star = 2 / d
    
    A = S_star * (beta_d / (gamma + d) + (beta_i * p1) / (Viral(0.114, 3.7594, current_temp) * (gamma + d)))
    B = S_star * (beta_d / (gamma + d + d2) + (beta_i * p2) / (Viral(a, b, current_temp) * (gamma + d + d2)))
    term1 = (1 - alpha1) * A + (1 - alpha2) * B
    term2 = (1 - alpha2) * A + (1 - alpha1) * B
    disc = term2**2 - 4 * (1 - alpha1) * (1 - alpha2) * A * B

    R0 = 0.5 * (term1 + np.sqrt(np.maximum(disc, 0)))  
    
    
    print(A, B, R0)
    import matplotlib.pyplot as plt 
    
    plt.figure( figsize = (10, 8))
    plt.plot(current_temp, A, color = 'b', label = '$R_0$ for LPAI' )
    plt.plot(current_temp, B, color = 'r', label = '$R_0$ for HPAI')
    plt.plot(current_temp, R0, color = 'black', alpha = 0.6, label = '$R_0$ for coexisting environment')
    plt.xlabel("Temperature (°C, sorted)")
    plt.ylabel("$R_0$")
    plt.title("Relationship between Temperature and $R_0$")
    plt.legend()
    plt.show()
    
    return R0
