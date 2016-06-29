import string
import math
import numpy as np
import matplotlib.pyplot as plt
    
data=open('../result/Ma=0.3/alpha=3/pressure_distribution.dat')

up_x=[]
up_p=[]

down_x=[]
down_p=[]

for line in data:
    if line=='\n':
        break
    (x,y,cp)=line.split()
    up_x.append(x)
    up_p.append(cp)

for line in data:
    (x,y,cp)=line.split()
    down_x.append(x)
    down_p.append(cp)

plt.scatter(up_x,up_p,label="UP")
plt.scatter(down_x,down_p,label="DOWN")
plt.axis([0, 1, 2, -2])
plt.xlabel("x_pos")
plt.ylabel("Cp")
plt.legend()
plt.show()

