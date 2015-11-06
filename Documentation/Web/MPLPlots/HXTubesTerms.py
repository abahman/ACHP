import matplotlib,pylab, numpy as np
from matplotlib.patches import FancyArrowPatch

def circle(x,y,r):
    """ This function makes lines for a circle given center and radius"""
    t=np.linspace(0,2*np.pi)
    xv=x+r*np.cos(t)
    yv=y+r*np.sin(t)
    return xv,yv

fig=pylab.figure()

## Actually make the set of tubes
nR=2
nC=1
offset = 0.5
r=0.25
for i in range(nR):
    for j in range(nC):
        if j%2==0:
            x=j
            y=i
        else:
            x=j
            y=i+offset
        xv,yv=circle(x,y,r)
        pylab.plot(xv,yv,'b')


##Dimension lines
pylab.plot(np.r_[-.25,-.25],np.r_[1,1.5],'k')
pylab.plot(np.r_[0.25,0.25],np.r_[1,1.5],'k')
pylab.text(0,1.5,'Tube\nDepth\n$T_d$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((-0.25,1.3),(0.25,1.3),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.plot(np.r_[-0.5,0],np.r_[0.75,0.75],'k')
pylab.plot(np.r_[-0.5,0],np.r_[0.25,0.25],'k')
pylab.text(-0.6,0.5,'Tube\nSpacing\n$b$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((-0.3,0.25),(-0.3,0.75),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.plot(np.r_[0.5,0],np.r_[0.25,0.25],'k')
pylab.plot(np.r_[0.5,0],np.r_[-0.25,-0.25],'k')
pylab.text(0.6,0.0,'Major\nDiameter\n$H_T$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((0.3,-0.25),(0.3,0.25),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

#Airflow arrow
pylab.gca().add_patch(FancyArrowPatch((0,0.5),(1.25,0.5),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(0.5,0.5,'Airflow\nDirection',ha='center',va='center')

pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.savefig('images/HXTubesTerms.pdf')
pylab.show()