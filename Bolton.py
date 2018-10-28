
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


C_to_K = 273.15
c_p_dry = 1005
c_v_dry = 1875
eps = .622
k_dry = .2854


# In[3]:


def sat_vapor_pressure(T):
    """vapor pressure for saturated air in mb as it varies with temperature"""
    e_s = 6.112 * np.exp((17.67 * T) / (T + 243.5))
    return e_s

def sat_vapor_temperature(e_s):
    """vapor temperature in Celsius as it varies with saturation vapor temperature"""
    e_t = ((243.5 * np.log(e_s) - 440.8) / (19.48 - np.log(e_s)))
    return e_t

def sat_mixing_ratio(p,T):
    """saturation mixing ratio in kg/kg"""
    w_s = eps * (sat_vapor_pressure(T)/(p - sat_vapor_pressure(T)))
    return w_s

# In[9]:


def mixing_ratio_line(p, w_s):
    """mixing ratio line that takes in pressure in mb sat_mixing_ratio kg/kg and returns a temperature in Celsius"""
    w_line = ((w_s * p) / (w_s + eps))
    return w_line


# In[10]:


def RH(T, p, w):
    """calculating RH in percent for a given temperature in Celsius and pressure in mb and mixing ratio kg/kg"""
    w_s = sat_mixing_ratio(p, T)
    RH = (w / w_s) * 100
    return RH


# In[11]:


def T_LCL(T, w, p, p_0=1000):
    """Given a temperature in Kelvin and RH in percent we will return a temperature at the LCL in Kelvin"""
    temp_lcl = (1 / (1 / ((T + C_to_K) - 55)) - (np.log(RH(T, p, w) / 100) / 2840)) + 55
    return temp_lcl
print(T_LCL(273.15, .0148836, 1000))

# In[12]:


def theta_dry(theta, p, p_0=1000.0):
    """Given a potential temperature theta in Kelvin and many possible pressures in mb return a temperature in Kelvin for potential dry air"""
    T_k = theta / (1000 / p)**k_dry
    return T_k


# In[20]:


def pseudoep_potential_T(T, w, p, p_0=1000.0):
    """given temperature in Celsius and pressure in mb and mixing ratio in both g/kg and kg/kg return theta_ep in Kelvin which is relative to moisture content of the air"""
    rh = RH(T, p, w)
    temp_lcl = T_LCL(T, w, p, p_0=1000)
    pseudoep_potential_T = (T + C_to_K) * (p_0 / p)**(0.2854 * (1 - (0.28 *  w))) * np.exp(((3.376 / temp_lcl) - 0.00254) * ((10**-3 * w) * (1 + (0.81 * w))))
    return pseudoep_potential_T
#print(pseudoep_potential_T(20, .0148836, 1000))

# In[21]:


def theta_ep_field(T, p, p_0=1000.0):
    """Similarly to the theta_e_field, we are writing out an equation that takes in any givnen temperature in Celsius and pressure in mb and we will return our moist adiabats for our Skew_T_Log_p diagram"""
    w_s = sat_mixing_ratio(p, T)
    moist_adiabat = pseudoep_potential_T(T, p, w_s, p_0=1000)
    return moist_adiabat
