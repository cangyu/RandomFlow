import string
import math
import numpy as np
import matplotlib.pyplot as plt

data=open('../results/Test1_1.txt')

for line in data:
    if line=="\n":
        break

data.readline()

index=[]
X=[]
Fx=[]
fitVal=[]

for line in data:
    (i,x,fx,fit)=line.split()
    index.append(int(i))
    X.append(float(x))
    Fx.append(math.log10(float(fx)))
    fitVal.append(float(fit))

plt.plot(index,X)
plt.xlabel("Iteration Step")
plt.ylabel("X")
plt.plot()
plt.savefig('../pic/T1-1_res_x.png')
plt.close()

plt.plot(index,Fx)
plt.xlabel("Iteration Step")
plt.ylabel("log(F(x))")
plt.plot()
plt.savefig('../pic/T1-1_res_fx.png')
plt.close()

