import string
import math
import numpy as np
import matplotlib.pyplot as plt

Mach=[0.3,0.8,1.2]
alpha=[0,3,6,9]

for Ma in Mach:
    for angle in alpha:

        up={}
        down={}

        input_path="../result/"+str(Ma)+"_"+str(angle)+"_"+"pressure_distribution.dat"
        output_path="../pic/"+str(Ma)+"_"+str(angle)+"_"+"Cp.png"

        data=open(input_path)

        for line in data:
            if line=='\n':
                break
            (x,y,cp)=line.split()
            up[float(x)]=cp

        for line in data:
            (x,y,cp)=line.split()
            down[float(x)]=cp

        up=sorted(up.items(),key=lambda d:d[0])
        down=sorted(down.items(),key=lambda d:d[0])

        up_x=[]
        up_p=[]

        down_x=[]
        down_p=[]

        for pos,cp in up:
            up_x.append(pos)
            up_p.append(cp)

        for pos,cp in down:
            down_x.append(pos)
            down_p.append(cp)

        plt.plot(up_x,up_p,label="UP")
        plt.plot(down_x,down_p,label="DOWN")
        plt.xlabel("X")
        plt.ylabel("Cp")
        plt.legend()
        plt.gca().invert_yaxis()
        plt.savefig(output_path)
        plt.close()
