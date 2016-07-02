import string
import math
import numpy as np
import matplotlib.pyplot as plt
import re

data=open('../results/Test2(single-thread).txt')

for line in data:
    if line=="\n":
        break

data.readline()

index=[]
p1=[]
p2=[]

for line in data:
    if line=="\n":
        continue

    t=line.split()
    if len(t)!=5:
        break
    
    if t[4]=='0':
        p1.append(t[3])
    elif t[4]=='1':
        p2.append(t[3])
    else:
        index.append(t[0]) 

plt.scatter(index,p1,label="Player1",color='m',marker='x')
plt.scatter(index,p2,label="Player2",color='r',marker='o')
plt.xlabel("Iteration Step")
plt.ylabel("F(x1,x2)")
plt.legend()
plt.plot()
plt.savefig('../pic/T2_res_f.png')
plt.close()


