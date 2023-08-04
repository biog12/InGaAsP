#--------------------------

# 已知上一炉 x,1-x,TmGa,TMIn
# 求下一炉,TMGa的使用量
# x/(1-x)=TMGa/TMIn

# # 已知上一炉 y,1-y,AsH3,PH3
# # 求下一炉,AsH3的使用量
#y/(1-y)=AsH3/PH3
#--------------------------

# a=旧x,c=新x，b=旧Ga
def f(a,b,c):
    x1=a
    x11=1-a
    x2=c
    x22=1-c
    return((x2/x22)/(x1/x11)*b)

x=f(0.256,7.9,0.26)
y=f(0.523,33.5,0.565)
print(x)
print(y)
