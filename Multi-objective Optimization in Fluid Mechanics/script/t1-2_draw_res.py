import string
import math
import numpy as np
import matplotlib.pyplot as plt

data=open('../results/Test1_2.txt')

for line in data:
    if line=="\n":
        break

data.readline()

index=[]
X1=[]
X2=[]
F=[]
fitVal=[]

for line in data:
    (i,x1,x2,f,fit)=line.split()
    index.append(int(i))
    X1.append(float(x1))
    X2.append(float(x2))
    F.append(math.log10(float(f)))
    fitVal.append(float(fit))

plt.plot(index,X1,label="X1")
plt.plot(index,X2,label="X2")
plt.xlabel("Iteration Step")
plt.ylabel("X")
plt.legend()
plt.plot()
plt.savefig('../pic/T1-2_res_x1_x2.png')
plt.close()

plt.plot(index,F)
plt.xlabel("Iteration Step")
plt.ylabel("log(F(x1,x2))")
plt.plot()
plt.savefig('../pic/T1-2_res_f.png')
plt.close()

