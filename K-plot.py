import numpy as np
import matplotlib.pyplot as plt

# y solid ratio
# x gas ratio
# k distribution factor
def function(x, k):
    return (k*x)/(1-x+k*x)

# 设置x的取值范围
x = np.linspace(0, 1, 100)

# 设置k的值
k = 0.04

# 计算对应于每个x值的y值
y = function(x, k)

# 绘制函数图像
plt.plot(x, y,color='red')
plt.xlabel('x')
plt.ylabel('y')
plt.title('Function: y / (1 - y) - (k * x) / (1 - x)')
plt.grid(True)
plt.show()