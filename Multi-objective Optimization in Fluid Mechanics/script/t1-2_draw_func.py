import numpy as np  
import matplotlib.pyplot as plt  
import mpl_toolkits.mplot3d
  
x,y=np.mgrid[-3:3:60j,-3:3:60j]  
z=(x-1.5)**4+(y+1.5)**2
 
ax=plt.subplot(111,projection='3d')  
ax.plot_surface(x,y,z,rstride=2,cstride=1,cmap=plt.cm.coolwarm,alpha=0.8)  
ax.set_xlabel('x1')  
ax.set_ylabel('x2')  
ax.set_zlabel('F(x1,x2)')  
  
plt.show()  