# -*- coding: utf-8 -*-
"""
Created on Mon Aug 28 22:11:25 2023

@author: biog1
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Jul 15 22:49:40 2023

@author: biog1
"""
#In(1-x-y)GaxAlyAs

# Constants from tab.1  JOQE 35/5 1999

a0_GaAs = 0.56533  
a0_InAs = 0.60584
a0_InP = 0.58688
a0_GaP = 0.54505
a0_AlAs = 0.566

c11_GaAs = 11.879
c11_InAs = 8.329
c11_InP = 10.110  
c11_GaP = 14.050
c11_AlAs = 12.500

c12_GaAs = 5.376
c12_InAs = 4.526
c12_InP = 5.610
c12_GaP = 6.203
c12_AlAs = 5.340

ac_GaAs = -7.17
ac_InAs = -5.08
ac_InP = -5.04
ac_GaP = -7.14  
ac_AlAs = -5.64

av_GaAs = 1.16
av_InAs = 1.00
av_InP = 1.27
av_GaP = 1.70
av_AlAs = 2.47

b_GaAs = -1.7
b_InAs = -1.8  
b_InP = -1.7
b_GaP = -1.8
b_AlAs = -1.5

g1_GaAs = 6.8
g1_InAs = 20.4
g1_InP = 4.95
g1_GaP = 4.05
g1_AlAs = 3.45

g2_GaAs = 1.9
g2_InAs = 8.3
g2_InP = 1.65
g2_GaP = 0.49
g2_AlAs = 0.68

g3_GaAs = 2.73
g3_InAs = 9.1
g3_InP = 2.35
g3_GaP = 1.25
g3_AlAs = 1.29

#----------------------------------------------------
from scipy.optimize import bisect

# PQ_eps_xy_f1
def PQ_eps_xy_f1(x, y):
    return (a0_InP - PQ_interp_a0(x,y)) / PQ_interp_a0(x,y)

# eps_z    
def eps_z(x, y):
    return (-2*PQ_interp_c12(x,y)/PQ_interp_c11(x,y) * PQ_eps_xy_f1(x,y))

# dEc_f1
def dEc_f1(x, y):
    return (PQ_interp_ac(x,y) * (2 * PQ_eps_xy_f1(x,y) + eps_z(x,y)))

# Peps_f1  
def Peps_f1(x, y):
    return (-2*PQ_interp_av(x,y) * (1 - PQ_interp_c12(x,y) / PQ_interp_c11(x,y)) * PQ_eps_xy_f1(x,y))

# Qeps_f1
def Qeps_f1(x, y):
    return (-PQ_interp_b(x,y)*(1+2*PQ_interp_c12(x,y)/PQ_interp_c11(x,y)) * PQ_eps_xy_f1(x,y))

# dEhh_f1
def dEhh_f1(x, y):
    return (-Peps_f1(x,y)-Qeps_f1(x,y))

# dElh_f1
def dElh_f1(x, y):
    return (-Peps_f1(x,y)+Qeps_f1(x,y))

# PQ_Eg_unstrained
def PQ_Eg_unstrained(x, y):
    return (0.36+2.093*y+0.629*x+0.577*y*y+0.436*x*x+1.013*x*y-2.0*x*y*(1-x-y)) #eV

# # PQ_Eg_unstrained_f1
# def PQ_Eg_unstrained_f1(x, y):
#     return (1.35+0.668*x-1.068*y+0.758*x*x+0.078*y*y-0.069*x*y-0.332*x*x*y+0.03*x*y*y) #eV

# PQ_interp_a0
def PQ_interp_a0(x, y):
    return (a0_GaAs*x + a0_AlAs*y + a0_InAs*(1-x-y))

# PQ_interp_c11
def PQ_interp_c11(x, y):
    return (c11_GaAs*x + c11_AlAs*y + c11_InAs*(1-x-y))

# PQ_interp_c12  
def PQ_interp_c12(x, y):
    return (c12_GaAs*x + c12_AlAs*y + c12_InAs*(1-x-y))

# PQ_interp_ac
def PQ_interp_ac(x, y):
    return (ac_GaAs*x + ac_AlAs*y + ac_InAs*(1-x-y))

# PQ_interp_av
def PQ_interp_av(x, y):
    return (av_GaAs*x + av_AlAs*y + av_InAs*(1-x-y))  

# PQ_interp_b
def PQ_interp_b(x, y):
    return (b_GaAs*x + b_AlAs*y + b_InAs*(1-x-y))

# PQ_Echh_f1
def PQ_Echh_f1(x, y):
    return (PQ_Eg_unstrained(x,y)+dEc_f1(x,y)-dEhh_f1(x,y))

# PQ_Eclh_f1  
def PQ_Eclh_f1(x, y):
    return (PQ_Eg_unstrained(x,y)+dEc_f1(x,y)-dElh_f1(x,y))

# PQ_Eg_XY
def PQ_Eg_XY(x, y):
    Eclh = PQ_Eclh_f1(x,y)
    Echh = PQ_Echh_f1(x,y)
    if Eclh <= Echh:
        return Eclh
    else:
        return Echh
        
# LatticeConst        
def LatticeConst(eps_xy):
    return (a0_InP/(eps_xy+1))
#*****************************************************
# PL2Eg_f1
#*****************************************************
def PL2Eg_f1(PL):
    h = 6.62606957e-34  # Planck constant
    c = 299792458   # speed of light
    q = 1.602176565e-19 # electron charge
    return (h*c/PL/q*1e9) # in eV

#*****************************************************
# Eg2PL_f1  
#***************************************************** 
def Eg2PL_f1(Eg):
    h = 6.62606957e-34  # Planck constant 
    c = 299792458   # speed of light
    q = 1.602176565e-19 # electron charge
    return (1e9*h*c/Eg/q) # PL in nm

#*****************************************************










# 全局变量  
EG2MATCH = None  
PQ_EPS_XY = None

# 计算晶格常数
def LATTICE_CONST(eps_xy):
  return a0_InP / (eps_xy + 1)

# 根据X和应变计算Y
def PQ_Y_XSTRAIN(x, eps_xy):
  a0 = LATTICE_CONST(eps_xy)
  nom = x * (a0_InP - a0_GaP) + a0 - a0_InP
  denom = x * (a0_GaAs - a0_GaP - a0_InAs + a0_InP) + a0_InAs - a0_InP
  return nom / denom

# 定义根函数
def FUNC4ROOT(x):
  y = PQ_Y_XSTRAIN(x, PQ_EPS_XY)
  return PQ_Eg_XY(x, y) - EG2MATCH
  
# 二分法终止条件
def ROOT_TERMINATION(min_x, max_x):
  return abs(max_x - min_x) <= 1e-6  

# 主函数
def PQ_X_EGSTRAIN(eg, eps_xy):

  global EG2MATCH, PQ_EPS_XY
  EG2MATCH = eg
  PQ_EPS_XY = eps_xy

  # 二分法寻根 
  result = bisect(FUNC4ROOT, 0, 1, xtol=1e-8)
  return result


#*****************************************************
# cal x,y from lamda eps + for compress -for tensil
#*****************************************************
def PQ_XY_Eg_strain(PL,eps):
    eps=-eps
    eg=1240.0/PL
    x=PQ_X_EGSTRAIN(eg,eps)
    y=PQ_Y_XSTRAIN(x,eps)
    return (x,y)


#*****************************************************
# cal  Eg eps from x,y,eps + for compress -for tensil
#*****************************************************
def PQ_lamda_eps_XY(x,y): 
    eg=PQ_Eg_XY(x,y)
    lamda=1240.0/eg
    eps=PQ_eps_xy_f1(x,y)
    return (lamda,-eps)
#*****************************************************




# Create a meshgrid of x and y values
import numpy as np
import matplotlib.pyplot as plt
import streamlit as st
from ternary_diagram import TernaryDiagram
# Define the x and y arrays
x = [0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
y = [0.0, 0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]

results={}


for x_val in x:
    for y_val in y:
         if (x_val+y_val)<=1.0:
            wavelengths, strains_xy = PQ_lamda_eps_XY(x_val, y_val)
            results[(x_val, y_val,1-x_val-y_val)] = (wavelengths, strains_xy)

#取坐标
keys = [key for key in results]
#转换为数组
xyz_array = np.array(keys)

x=xyz_array[:,0]
y=xyz_array[:,1]
z=xyz_array[:,2]

#取值
values = [results[key] for key in results]

v_array=np.array(values)

PL=v_array[:,0]
Strain=v_array[:,1]
td = TernaryDiagram(["GaAS", "AlAs", "InAs"])
Level_PL =[1100,1300]

td.contour(xyz_array,PL,z_min=1000,z_max=1500)


plt.show()
