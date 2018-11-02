
# coding: utf-8

# In[1]:


import numpy as np


# In[2]:


C_to_K = 273.15  # J/kg*K
c_p_dry = 1005  #J/kg*K
c_v_dry = 1875   #J/kg*K
eps = .622   #dimensionless
k_dry = .2854   #dimensionless


# In[3]:


def sat_vapor_pressure(T):
    """vapor pressure for saturated air in mb as it varies with temperature in Celsius"""
    e_s = 6.112 * np.exp((17.67 * (T - C_to_K)) / ((T - C_to_K) + 243.5))
    return e_s
# print(sat_vapor_pressure(20))

def sat_vapor_temperature(e_s):
    """vapor temperature in Celsius as it varies with saturation vapor temperature"""
    e_t = (243.5 * np.log(e_s) - 440.8) / (19.48 - (np.log(e_s)))
    return e_t

def sat_mixing_ratio(T, p):
    """saturation mixing ratio in kg/kg where pressure is in Pascals"""
    e_s = sat_vapor_pressure(T - C_to_K)
    w_s = eps * (e_s / ((p) - e_s))
    return w_s
# print(sat_mixing_ratio(p, T))
# In[9]:


def mixing_ratio_line(p, w_s):
    """mixing ratio line that takes in pressure in mb sat_mixing_ratio kg/kg and returns a temperature in Celsius"""
    w_l = ((w_s * (p)) / (eps + w_s))
    return w_l


# In[10]:


def RH(T, p, w):
    """calculating RH in percent for a given temperature in Celsius and pressure in mb and mixing ratio kg/kg"""
    w_s = sat_mixing_ratio(T - C_to_K, (p / 100))
    rh = ((w) / w_s) * 100
    return rh
#print(RH(20, 1000, .0111))

# In[11]:


def T_LCL(T, p, w, p_0=1000.0):
    """Given a temperature in Kelvin and RH in percent we will return a temperature at the LCL in Kelvin"""
    t_lcl = (1 / (1 / (T - 55)) - ((np.log(RH(T - C_to_K, p / 100, w) / 100) / 2840))) + 55
    return t_lcl

# In[12]:


def theta_dry(theta, p, p_0=1000.0):
    """Given a potential temperature theta in Kelvin and many possible pressures in mb return a temperature in Kelvin for potential dry air"""
    T_k = theta * (1000.0 / (p / 100))**k_dry
    return T_k


# In[20]:


def pseudoep_potential_T(T, p, w, p_0=1000.0):
    """given temperature in Celsius and pressure in mb and mixing ratio in both g/kg and kg/kg return theta_ep in Kelvin which is relative to moisture content of the air"""
    rh = RH(T, p / 100, w)
    t_lcl = T_LCL(T, p / 100, w, p_0=1000.0)
    theta_ep = (T) * (1000.0 / (p / 100))**(0.2854 * (1 - (0.28 *  w))) * np.exp((3.376 / t_lcl) - 0.00254) * ((w / 1000) * (1 + (0.81 * w)))
    return theta_ep
#print(pseudoep_potential_T(20, .0148836, 1000))

# In[21]:


def theta_ep_field(T, p, p_0=1000.0):
    """Similarly to the theta_e_field, we are writing out an equation that takes in any givnen temperature in Celsius and pressure in mb and we will return our moist adiabats for our Skew_T_Log_p diagram"""
    w_s = sat_mixing_ratio(T, p / 100)
    moist_adiabats = pseudoep_potential_T(T, p / 100, w_s, p_0=100000.0)
    return moist_adiabats
