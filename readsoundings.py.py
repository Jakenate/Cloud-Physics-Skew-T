
# coding: utf-8

# In[1]:


import numpy as np


# In[3]:


def parse_SPC("MAF-sounding data.txt", skip_rows=6):
    dtype = [('p', float), ('z', float), ('T', float), ('Td', float), ('wind_dir', float), ('wind_spd', float)]
    data = np.genfromtxt("MAF-sounding data.txt", dtype=dtype, skip_header=skip_rows, delimiter=',')
    return data
