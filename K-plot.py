import numpy as np
import matplotlib.pyplot as plt

# y solid ratio
# x gas ratio
# k distribution factor
def function(x, k):
    return (k*x)/(1-x+k*x)

def wl(y):
    Eg=1.35-0.775*y+0.149*y*y
    wl=1240.0/Eg
    return(wl)

# 设置x的取值范围
x = np.linspace(0, 1, 100)

# 设置k的值
k = 0.026

# 计算对应于每个x值的y值
y = function(x, k)
y2=wl(y)

# 绘制函数图像
fig,ax1=plt.subplots()
ax2=ax1.twinx()
##
ax1.plot(x, y,color='red')
ax1.set_xlabel('x-Gas phase')
ax1.set_ylabel('y-Solid phase')
##
ax2.plot(y, y2,color='black')
ax2.set_ylabel('wavelength(nm)')

plt.title('Gas vs Solid Phase and Wavelength')
plt.grid(True)
plt.show()