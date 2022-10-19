import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation 
import time

def golden_section_update():
    
    global a,b,c,d,q

    # compute interior points
    c = b - (b-a)/phi 
    d = a + (b-a)/phi 
    
    if f(c) < f(d):
        
        # point to remove
        q = b
        
        # update
        b = d
        
    else:
        
        # point to remove
        q = a
        
        # update
        a = c


# minimize te function from section 7
f = lambda x: np.exp(0.5*x-1)*(x+1)**2

# initial bracket
a = -8 
b = 1

# stopping tolerence (bracet size)
tolerence = 1e-8

# maximum number of iterations
max_iterations = 30

# golden ratio
phi = (1+np.sqrt(5))/2

# create plot 
fig,ax = plt.subplots(1)
ax.set_xlabel('x',fontsize = 16); ax.set_ylabel('f(x)',fontsize = 16)

x_ = np.linspace(a,b,10000)

ax.plot(x_,f(x_),linewidth = 2)

def animate(i):

    title = 'Iteration {:02d}'.format(i)
    
    # golden section update
    golden_section_update()

    # set picture
    ax.set_title(title)
    
    # plot interior points
    ax.plot(c,f(c),'r.',markersize = 30)
    ax.plot(d,f(d),'g.',markersize = 30)
            
    # update plot with removed point
    ax.plot(q,f(q),'.',markersize=30,color = [0.5,0.5,0.5])
    
    return ax

anim = FuncAnimation(fig,animate,frames = list(range(max_iterations)),interval=500,blit=False,repeat=False)
anim.save("golden_section.gif", dpi=120, writer="imagemagick")

x = (a+b)/2

print(x)
            
