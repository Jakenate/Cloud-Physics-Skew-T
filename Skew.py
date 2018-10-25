
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import FuncFormatter, Formatter
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

# In[2]:


C_to_K = 273.15


# In[3]:


skew_slope = 40


# In[4]:


def x_from_Tp(T, p):
    """Transform x coordinate to temperature in degrees Celsius and pressure in mb"""
    x = T - skew_slope * (np.log(p))
    return x

def y_from_p(p):
    """Transform y coordinate to pressure in mb"""
    y = -(np.log(p))
    return y

def T_from_xp(x, p):
    """transform temperature to get back x coordinate and pressure in mb"""
    T = x + skew_slope*(np.log(p))
    return T

def p_from_y(y):
    """transform pressure in mb to get back y coordinate"""
    p = (np.exp(-y))


# In[5]:


def to_thermo(x, y):
    """transform (x, y) coordinates to T in degrees Celsius and p in mb."""
    p = p_from_y(y)
    T_C = T_from_xp(x, p) - C_to_K
    return p, T_C

def from_thermo(p, T_C):
    """transform T_C (in degrees Celsius) and p (in mb) to (x, y)"""
    y = y_from_p(p)
    x = x_from_Tp(T_C + C_to_K, p)
    return x, y


# In[6]:


# values along the botttom and left edges
p_bottom = 1050.0
p_top = 150
T_min = -40 + C_to_K
T_max = 50 + C_to_K


# In[7]:


x_min = np.min(x_from_Tp(T_min, p_bottom))
x_max = np.max(x_from_Tp(T_max, p_top))
y_min = np.min(y_from_p(p_bottom))
y_max = np.max(y_from_p(p_top))


# In[8]:


print(x_min)
print(x_max)
print(y_min)
print(y_max)


# In[9]:


p_levels = np.arange(1000, 100 - 50, -50)


# In[10]:


print(p_levels)


# In[11]:


T_C_levels = np.arange(-80, 40 + 10, 10)


# In[12]:


print(T_C_levels)


# In[13]:


T_levels = T_C_levels + C_to_K


# In[14]:


theta_levels = np.arange(-40, 100 + 10, 10) + C_to_K


# In[15]:


print(theta_levels)


# In[16]:


theta_ep_levels = theta_levels.copy()


# In[17]:


mixing_ratios = np.asarray([.4, 1, 2, 3, 5, 8, 12, 16, 20])


# In[ ]:


import Bolton

p_all = np.arange(p_bottom, p_top - 1, -1)

y_p_levels = y_from_p(p_levels)

y_all_p = y_from_p(p_all)

x_T_levels = [x_from_Tp(Ti, p_all) for Ti in T_levels]

x_thetas = [x_from_Tp(Bolton.theta_dry(theta_i, p_all), p_all) for theta_i in theta_levels]

x_mixing_ratios = [x_from_Tp(Bolton.mixing_ratio_line(p_all, mixing_ratio_i) + C_to_K, p_all) for mixing_ratio_i in mixing_ratios]

mesh_T, mesh_p = np.meshgrid(np.arange(-60.0, T_levels.max() - C_to_K + 0.1, 0.1), p_all)
theta_ep_mesh = Bolton.theta_ep_field(mesh_T, mesh_p)

def theta_e(T, p, p_0=1000):
     R_d = 287.058
     alpha = 3.139 * 10**6
     c_p = 1005
     c_l = 4218
     c_wd = c_p + (w_s * c_l)
     L_v = alpha - (c_l - c_p) * T
     w_s = sat_mixing_ratio(p, T)
     theta_e = T * (p_0 / p)**(R_d / c_wd) * (np.exp((L_v * w_s) / (c_wd * T)))
     return theta_e

def theta_e_field(T, p, p_0=1000.0):
    t_e = theta_e_field(T, p, p=1000.0)
    theta_e_field = t_e
    return theta_e_field

fig = plt.figure()
ax = Subplot(fig, 1, 1, 1)
fig.add_subplot(ax)

for yi in y_p_levels:
    ax.plot((x_min, x_max), (yi, yi), color=(1.0, 0.8, 0.8))

for x_T in x_T_levels:
    ax.plot(x_T, y_all_p, color=(1.0, 0.5, 0.5))

for x_theta in x_thetas:
    ax.plot(x_theta, y_all_p, color=(1.0, 0.7, 0.7))

for x_mixing_ratio in x_mixing_ratios:
    good = p_all >= 600 #restrict mixing ratio lines to below 600 mb
    ax.plot(x_mixing_ratio[good], y_all_p[good], color=(0.8, 0.8, 0.6))

n_moist = len(theta_ep_levels)
moist_colors = ((0.6, 0.9, 0.7),)*n_moist
ax.contour(x_from_Tp(mesh_T + C_to_K, mesh_p), y_from_p(mesh_p), theta_ep_mesh, theta_ep_levels, colors=moist_colors)



ax.axis((x_min, x_max, y_min, y_max))

plt.show()
