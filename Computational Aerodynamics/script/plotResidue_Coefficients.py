import string
import math
import numpy as np
import matplotlib.pyplot as plt
    
data=open('../result/Ma=0.3/alpha=0/convergence_history.dat')

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
plt.savefig('../pic/Residual_0.3_0.png')
plt.close()

plt.plot(index,CL,label="CL")
plt.plot(index,CD,label="CD")
plt.xlabel("Iteration Step")
plt.ylabel("Aerodynamic coefficients")
plt.legend()
plt.savefig('../pic/Coefficient_0.3_0.png')
plt.close()

