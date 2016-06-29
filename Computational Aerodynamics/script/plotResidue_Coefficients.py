import string
import math
import numpy as np
import matplotlib.pyplot as plt

Mach=[0.3,0.8,1.2]
alpha=[0,3,6,9]

for Ma in Mach:
    for angle in alpha:

        input_path="../result/"+str(Ma)+"_"+str(angle)+"_"+"convergence_history.dat"
        err_output_path="../pic/"+str(Ma)+"_"+str(angle)+"_"+"residual.png"
        coeff_output_path="../pic/"+str(Ma)+"_"+str(angle)+"_"+"cl_cd.png"
        
        data=open(input_path)

        index=[]
        err_rho=[]
        err_vel=[]
        err_p=[]

        CL=[]
        CD=[]

        for record in data:
            (i,er,ev,ep,cl,cd)=record.split()
            er=math.log10(float(er))
            ev=math.log10(float(ev))
            ep=math.log10(float(ep))
            index.append(i)
            err_rho.append(er)
            err_vel.append(ev)
            err_p.append(ep)
            CL.append(cl)
            CD.append(cd)
    
        plt.plot(index,err_rho,label="err_rho")
        plt.plot(index,err_vel,label="err_vel")
        plt.plot(index,err_p,label="err_p")
        plt.xlabel("Iteration Step")
        plt.ylabel("log(err)")
        plt.legend()
        plt.savefig(err_output_path)
        plt.close()

        plt.plot(index,CL,label="CL")
        plt.plot(index,CD,label="CD")
        plt.xlabel("Iteration Step")
        plt.ylabel("Aerodynamic coefficients")
        plt.legend()
        plt.savefig(coeff_output_path)
        plt.close()

