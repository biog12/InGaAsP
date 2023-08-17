#--------------------------

# 已知上一炉 x,1-x,TmGa,TMIn
# 求下一炉,TMGa的使用量
# x/(1-x)=TMGa/TMIn

# # 已知上一炉 y,1-y,AsH3,PH3
# # 求下一炉,AsH3的使用量
#y/(1-y)=AsH3/PH3
#--------------------------


import streamlit as st

# Function to calculate new TMGa
def calculate_new_tmgas(tmgas, ga_composition):
    old_x = ga_composition
    old_1_x = 1 - ga_composition
    new_x = st.number_input('Enter new Ga composition:', min_value=0.1, max_value=20.0, value=ga_composition)
    new_1_x = 1 - new_x
    new_tmgas = (new_x / new_1_x) / (old_x / old_1_x) * tmgas
    return new_tmgas

# Function to calculate new AsH3
def calculate_new_ash3(ash3, as_composition):
    old_y = as_composition
    old_1_y = 1 - as_composition
    new_y = st.number_input('Enter new As composition:', min_value=0.1, max_value=100.0, value=as_composition)
    new_1_y = 1 - new_y
    new_ash3 = (new_y / new_1_y) / (old_y / old_1_y) * ash3
    return new_ash3

# Streamlit app
st.title('Gas Composition Calculator')

st.header('TMGa Calculation')
tmgas_input = st.number_input('Enter original TMGa:', value=3.0)
ga_composition_input = st.number_input('Enter original Ga composition:', min_value=0.1, max_value=20.0, value=3.0)
new_tmgas_output = calculate_new_tmgas(tmgas_input, ga_composition_input)
st.write('New TMGa:', new_tmgas_output)

st.header('AsH3 Calculation')
ash3_input = st.number_input('Enter original AsH3:', value=10.0)
as_composition_input = st.number_input('Enter original As composition:', min_value=0.1, max_value=100.0, value=10.0)
new_ash3_output = calculate_new_ash3(ash3_input, as_composition_input)
st.write('New AsH3:', new_ash3_output)