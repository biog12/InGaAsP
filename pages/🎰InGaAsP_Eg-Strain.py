# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 21:43:05 2023

@author: biog1
"""
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
    return (1.35+0.642*x-1.101*y+0.758*x*x+0.101*y*y-0.159*x*y-0.28*x*x*y+0.109*x*y*y) #eV

# PQ_Eg_unstrained_f1
def PQ_Eg_unstrained_f1(x, y):
    return (1.35+0.668*x-1.068*y+0.758*x*x+0.078*y*y-0.069*x*y-0.332*x*x*y+0.03*x*y*y) #eV

# PQ_interp_a0
def PQ_interp_a0(x, y):
    return (a0_GaAs*x*y + a0_GaP*x*(1-y) + a0_InAs*(1-x)*y + a0_InP*(1-x)*(1-y))

# PQ_interp_c11
def PQ_interp_c11(x, y):
    return (c11_GaAs*x*y + c11_GaP*x*(1-y) + c11_InAs*(1-x)*y + c11_InP*(1-x)*(1-y))

# PQ_interp_c12  
def PQ_interp_c12(x, y):
    return (c12_GaAs*x*y + c12_GaP*x*(1-y) + c12_InAs*(1-x)*y + c12_InP*(1-x)*(1-y))

# PQ_interp_ac
def PQ_interp_ac(x, y):
    return (ac_GaAs*x*y + ac_GaP*x*(1-y) + ac_InAs*(1-x)*y + ac_InP*(1-x)*(1-y))

# PQ_interp_av
def PQ_interp_av(x, y):
    return (av_GaAs*x*y + av_GaP*x*(1-y) + av_InAs*(1-x)*y + av_InP*(1-x)*(1-y))  

# PQ_interp_b
def PQ_interp_b(x, y):
    return (b_GaAs*x*y + b_GaP*x*(1-y) + b_InAs*(1-x)*y + b_InP*(1-x)*(1-y))

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
  return abs(max_x - min_x) <= 1e-8  

# 主函数
def PQ_X_EGSTRAIN(eg, eps_xy):

  global EG2MATCH, PQ_EPS_XY
  EG2MATCH = eg
  PQ_EPS_XY = eps_xy

  # 二分法寻根 
  result = bisect(FUNC4ROOT, 0, 1, xtol=1e-6)
  return result


#*****************************************************
# cal x,y from lamda eps + for compress -for tensil,so ep=-eps
#*****************************************************
def PQ_XY_Eg_strain(PL,eps):
    ep=-eps
    eg=1240.0/PL
    x=PQ_X_EGSTRAIN(eg,ep)
    y=PQ_Y_XSTRAIN(x,ep)
    return (x,y)


#*****************************************************
# cal  Eg eps from x,y + for compress -for tensil,so ep=-eps
#*****************************************************
def PQ_lamda_eps_XY(x,y): 
    eg=PQ_Eg_XY(x,y)
    lamda=1240.0/eg
    eps=PQ_eps_xy_f1(x,y)
    return (lamda,-eps)
#*****************************************************

#######Ga(x)In(1-x)As(y)P(1-y) Bandgap Calculator#############################################################
import streamlit as st

st.title('Ga(x)In(1-x)As(y)P(1-y) Bandgap Calculator')

# Input widgets
pl = st.number_input('PL (nm)', value=1300)

# Input Epsilion for ppm
eps = st.number_input('Epsilon(ppm)', value=0.0, step=0.000001)/1000000

x = st.number_input('x', value=0.2, step=0.001, format='%.3f')
y = st.number_input('y', value=0.5, step=0.001, format='%.3f')

# Functions

# Create a column layout for the buttons
col1, col2 = st.columns(2)

with col1:
    if st.button('Calculate x,y'):
        x, y = PQ_XY_Eg_strain(pl, eps)
        st.write(f'x: {x:.3f}')
        st.write(f'y: {y:.3f}')

with col2:
    if st.button('Calculate PL,eps'):
        lamda, eps = PQ_lamda_eps_XY(x, y)
        eps=eps*1000000 #for ppm
        st.write(f'Wavelength: {lamda:.1f} nm')
        st.write(f'Epsilon: {eps:.1f} ppm')


