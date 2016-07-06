import string
import math
import numpy as np
import matplotlib.pyplot as plt

Re=[20,40,100,200]
results=[]
aerodynamic=[]

for re in Re:
    input_dir="../result/Re="+str(re)+"/aerodynamics.dat"
    output_dir="../pic/Aerodynamics_"+str(re)+".png"
    results.append(input_dir)
    aerodynamic.append(output_dir)

results.append("../result/Contrast/aerodynamics.dat")
aerodynamic.append("../pic/Aerodynamics_Contrast.png")

index=[]
Drag=[]
Lift=[]
cnt=0

for file in results:
    index.clear()
    Drag.clear()
    Lift.clear()

    data=open(file)
    for line in data:
        (i,cd,cl)=line.split()
        index.append(i)
        Drag.append(cd)
        Lift.append(cl)

    plt.plot(index,Drag,label="CD")
    plt.plot(index,Lift,label="CL")
    plt.legend()
    plt.xlabel("Iteration Step")
    plt.ylabel("Aerodynamic Coefficients")
    plt.savefig(aerodynamic[cnt])
    plt.close()
    data.close()
    cnt+=1

