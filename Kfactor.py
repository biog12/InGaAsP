# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 20:55:42 2023
InxGa1-xAS1-yPy

y/(1-y)=K(fAs/fP)

DOI:10.1016/0022-0248(96)00295-3

The P/As distribution coefficient K for growth of InGaAsP from 
phosphine and arsine is strongly dependent on temperature and 
the In/Ga content of the solid, and weekly dependent on other 
processing conditions such as total pressure
@author: biog12
"""


import numpy as np
from scipy.interpolate import interp2d
import streamlit as st

# Define the x and y arrays
x = [0, 0.25, 0.5, 0.75, 1]
y = [550, 600, 650, 700]

# Define the data array
data = [[0.0083, 0.02, 0.044, 0.088],
        [0.0073, 0.018, 0.038, 0.077],
        [0.0062, 0.015, 0.033, 0.066],
        [0.0052, 0.012, 0.027, 0.055],
        [0.0042, 0.01, 0.022, 0.044]]

# Create a meshgrid of x and y values
X, Y = np.meshgrid(x, y)

# Transpose the data array to match the X and Y arrays
Z = np.transpose(data)

# Create a 2D interpolation function
f = interp2d(x, y, Z)

# Streamlit app
st.title('P/As分配系数calculation App')

# Input temp and x_val
temp = st.number_input('输入温度 Temp:', value=645)
x_val = st.number_input('输入In_组分:', value=0.92)

# Calculate result using interpolation function
result = f(x_val, temp)

# Display result
# Calculate result using interpolation function
if st.button('Calculate'):
    result = f(x_val, temp)
    # Display result
    st.write('Result:', result)