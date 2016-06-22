import string
import math
import numpy as np
import matplotlib.pyplot as plt

Re=[20,40,100,200]
results=[]
residue=[]

for re in Re:
    input_dir="../result/Re="+str(re)+"/residue.dat"
    output_dir="../pic/Residue_"+str(re)+".png"
    results.append(input_dir)
    residue.append(output_dir)

results.append("../result/Contrast/residue.dat")
residue.append("../pic/Residue_Contrast.png")

index=[]
err_vel=[]
err_rho=[]
cnt=0

for file in results:
    index.clear()
    err_vel.clear()
    err_rho.clear()

    data=open(file)
    for line in data:
        (i,ev,er)=line.split()
        ev=math.log10(float(ev))
        er=math.log10(float(er))
        index.append(i)
        err_vel.append(ev)
        err_rho.append(er)

    plt.plot(index,err_vel)
    plt.plot(index,err_rho)
    plt.savefig(residue[cnt])
    plt.close()
    data.close()
    cnt+=1

