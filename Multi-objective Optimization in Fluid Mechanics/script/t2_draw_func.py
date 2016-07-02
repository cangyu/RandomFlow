import numpy as np  
import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
  
x,y=np.mgrid[-4:4:40j,-4:4:40j]  
z=(x*y-y)**2+(x*y-2)**4
 
ax=plt.subplot(111,projection='3d')  
ax.plot_surface(x,y,z,rstride=2,cstride=1,cmap=plt.cm.coolwarm,alpha=0.8)  
ax.set_xlabel('x1')  
ax.set_ylabel('x2')  
ax.set_zlabel('F(x1,x2)')  
  
plt.show()  
