
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
    e_s = 6.112 * np.exp((17.67 * T) / (T + 243.5))
    return e_s

def sat_vapor_temperature(e_s):
    e_t = ((243.5 * np.log(e_s) - 440.8) / (19.48 - np.log(e_s)))
    return e_t

def sat_mixing_ratio(p,T):
    w_s = eps * (e_s/(p - e_s))
    return w


# In[9]:


def mixing_ratio_line(p, w_s):
    w_line = ((w_s * p) / (w_s + eps))
    return w_line


# In[10]:


def RH(T, p, w):
    RH = (w / w_s) * 100
    return RH


# In[11]:


def T_LCL(T, RH):
    T_LCL = 1 / ((1 / (C_to_K - 55)) - ((np.log(RH / 100))/2840)) +55
    return T_LCL


# In[12]:


def theta_dry(theta, p, p_0 = 1000.0):
    theta_dry = (T + C_to_K) * (1000 / p)**k_dry
    return theta_dry


# In[20]:


def pseudoep_potential_T(T, theta_dry, w, p_0=1000.0):
    pseudoep_potential_T = (T + C_to_K) * (p_0 / p)**0.2854(1 - (0.28 * 10**-3 * w)) * np.exp(((3.376 / T_LCL) - 0.00254) * (w * (1 + (0.81 * (10)**-3 * w))))
    return pseudoep_potential_T


# In[21]:


def theta_ep_field(T, p, p_0=1000.0):
    w_s = sat_mixing_ratio(p, T)
    moist_adiabat = pseudoep_potential_T(T, p, w_s, p_0=p_0)
    return moist_adiabat

