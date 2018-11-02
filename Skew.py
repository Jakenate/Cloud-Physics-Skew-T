
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axisartist import Subplot
from matplotlib.ticker import FuncFormatter, Formatter
from mpl_toolkits.axisartist.grid_helper_curvelinear import GridHelperCurveLinear

# In[2]:


C_to_K = 273.15  #convert temperature in Celsius to Kelvin


# In[3]:


skew_slope = 40  #value for alpha in given equations


# In[4]:


def y_from_p(p):
    """Transform y coordinate to pressure in mb"""
    y = -(np.log(p))
    return y

def x_from_Tp(T, p):
    """Transform x coordinate to temperature in degrees Celsius and pressure in mb"""
    x = (T) - skew_slope * (np.log(p))
    return x

def p_from_y(y):
    """transform pressure in mb to get back y coordinate"""
    p = np.exp(-y)
    return p

def T_from_xp(x, p):
    """transform temperature to get back x coordinate and pressure in mb"""
    T = x - (skew_slope * y_from_p(p))
    return T


# In[5]:


def to_thermo(x, y):
    """transform (x, y) coordinates to T in degrees Celsius and p in mb."""
    p = p_from_y(y)
    T = T_from_xp(x, p) - C_to_K
    return T, p
    #print(T)

def from_thermo(p, T):
    """transform T_C (in degrees Celsius) and p (in mb) to (x, y)"""
    y = y_from_p(p)
    x = x_from_Tp(T + C_to_K, p)
    return x, y


# In[6]:


# values along the botttom and left edges
p_bottom = 105000.0  #in Pascals
p_top = 15000.0  #in Pascals
T_min = -40 + C_to_K  # in Kelvin
T_max = 50 + C_to_K  #in Kelvin


# In[7]:


x_min = np.min(x_from_Tp(T_min, p_bottom))
x_max = np.max(x_from_Tp(T_max, p_top))
y_min = np.min(y_from_p(p_top))
y_max = np.max(y_from_p(p_bottom))


# In[8]:


# print(x_min)
# print(x_max)
# print(y_min)
# print(y_max)


# In[9]:


p_levels = np.arange(100000, 10000 - 5000, -5000)  #in Pascals


# In[10]:


#print(p_levels)


# In[11]:


T_C_levels = np.arange(-80, 40 + 10, 10)  # in Kelvin


# In[12]:


#print(T_C_levels)


# In[13]:


T_levels = T_C_levels + C_to_K


# In[14]:


theta_levels = np.arange(-40, 100 + 10, 10) + C_to_K


# In[15]:


#print(theta_levels)


# In[16]:


theta_ep_levels = theta_levels.copy()


# In[17]:


mixing_ratios = np.asarray([.4, 1, 2, 3, 5, 8, 12, 16, 20]) / 1000.0


# In[ ]:


import Bolton

p_all = np.arange(p_bottom, p_top - 10000, -10000)
#print(p_all)

y_p_levels = y_from_p(p_levels)
#print(y_p_levels)

y_all_p = y_from_p(p_all)
#print(y_all_p)

x_T_levels = [x_from_Tp(Ti, p_all) for Ti in T_levels]

x_thetas = [x_from_Tp(Bolton.theta_dry(theta_i, p_all), p_all) for theta_i in theta_levels]

x_mixing_ratios = [x_from_Tp(Bolton.mixing_ratio_line(p_all, mixing_ratio_i) + C_to_K, p_all) for mixing_ratio_i in mixing_ratios]

mesh_T, mesh_p = np.meshgrid(np.arange(-60.0, T_levels.max() - C_to_K + 0.1, 0.1), p_all)
theta_ep_mesh = Bolton.theta_ep_field(mesh_T, mesh_p)

def ep_potential_T(T, p, p_0=1000.0):
    """equivalent potential temperature in Kelvin as it varies with temperature and pressure"""
    R_d = 287.058 #J/kg*K
    alpha = 3.139 * 10**6 #J/kg
    c_p = 1005 #J/kg*K
    c_l = 4218 #J/kg*K
    w_s = sat_mixing_ratio(T, p)
    c_wd = c_p + (w_s * c_l)
    L_v = alpha + (c_l - c_p) * T
    theta_e = T * (p_0 / p)**(R_d / c_wd) * (np.exp((L_v * w_s) / (c_wd * T)))
    return theta_e

def theta_e_field(T, p, p_0=1000.0):
    """create a theta_e field that varies with temperature and pressure"""
    t_e = theta_e_field(T, p, p=1000.0)
    theta_e_field = t_e
    return theta_e_field

skew_grid_helper = GridHelperCurveLinear((from_thermo, to_thermo))

fig = plt.figure()
ax = Subplot(fig, 1, 1, 1, grid_helper=skew_grid_helper)
fig.add_subplot(ax)

def format_coord(x, y):
    """format ticks with physical values"""
    T, p = to_thermo(x, y)
    return "{0: 5.1f} C {1: 5.1f} mb".format(float(T), float(p))

ax.format_coord = format_coord

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

#code for theta_e

ax.axis((x_min, x_max, y_min, y_max))

plt.show()
