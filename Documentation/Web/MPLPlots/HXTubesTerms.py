import matplotlib,pylab, numpy as np
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt

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
        #pylab.plot(xv,yv,'b')


verts1 = [
    (0.24, 0.05),  # P0
    (0.3, 0.05), # P1
    (0.3, -0.05), # P2
    (0.24, -0.05), # P3
    (-0.24, -0.05), # P4
    (-0.3, -0.05), # P5
    (-0.3, 0.05), # P6
    (-0.24, 0.05), ## P7
    (0.24, 0.05),  # to close the loop , please ignored because it is the same as the first point
    ]

verts2 = [
    (0.22, 0.035),  # P0
    (0.28, 0.035), # P1
    (0.28, -0.035), # P2
    (0.22, -0.035), # P3
    (-0.22, -0.035), # P4
    (-0.28, -0.035), # P5
    (-0.28, 0.035), # P6
    (-0.22, 0.035), ## P7
    (0.22, 0.035),  # to close the loop , please ignored because it is the same as the first point
    ]

verts3 = [
    (0.24, 1.05),  # P0
    (0.3, 1.05), # P1
    (0.3, 0.95), # P2
    (0.24, 0.95), # P3
    (-0.24, 0.95), # P4
    (-0.3, 0.95), # P5
    (-0.3, 1.05), # P6
    (-0.24, 1.05), ## P7
    (0.24, 1.05),  # to close the loop , please ignored because it is the same as the first point
    ]

verts4 = [
    (0.22, 1.035),  # P0
    (0.28, 1.035), # P1
    (0.28, 0.965), # P2
    (0.22, 0.965), # P3
    (-0.22, 0.965), # P4
    (-0.28, 0.965), # P5
    (-0.28, 1.035), # P6
    (-0.22, 1.035), ## P7
    (0.22, 1.035),  # to close the loop , please ignored because it is the same as the first point
    ]

codes = [Path.MOVETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         ]

path1 = Path(verts1, codes)
path2 = Path(verts2, codes)
path3 = Path(verts3, codes)
path4 = Path(verts4, codes)
ax = fig.add_subplot(111)
patch = patches.PathPatch(path1, facecolor='none',edgecolor='blue', lw=0.8)
patch2 = patches.PathPatch(path2, facecolor='none',edgecolor='blue', lw=0.8)
patch3 = patches.PathPatch(path3, facecolor='none',edgecolor='blue', lw=0.8)
patch4 = patches.PathPatch(path4, facecolor='none',edgecolor='blue', lw=0.8)

ax.add_patch(patch)
ax.add_patch(patch2)
ax.add_patch(patch3)
ax.add_patch(patch4)

#Uncomment the next two line to see the point's path used
#xs, ys = zip(*verts1)
#pylab.plot(xs, ys, 'x--', lw=0.8, color='black', ms=10)

##Dimension lines
pylab.plot(np.r_[-.3,-.3],np.r_[1,1.5],'k')
pylab.plot(np.r_[0.3,0.3],np.r_[1,1.5],'k')
pylab.text(0,1.4,'Tube\nDepth\n$T_d$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((-0.3,1.25),(0.3,1.25),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.plot(np.r_[-0.7,-0.3],np.r_[0.95,0.95],'k')
pylab.plot(np.r_[-0.7,-0.3],np.r_[0.05,0.05],'k')
pylab.text(-0.75,0.5,'Tube\nSpacing\n$b$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((-0.5,0.05),(-0.5,0.95),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.plot(np.r_[0.5,0.3],np.r_[0.05,0.05],'k')
pylab.plot(np.r_[0.5,0.3],np.r_[-0.05,-0.05],'k')
pylab.text(0.7,0.0,'Major\nDiameter\n$H_t$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((0.4,0.05),(0.4,0.2),arrowstyle='<|-',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.gca().add_patch(FancyArrowPatch((0.4,-0.2),(0.4,-0.05),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

#Airflow arrow
pylab.gca().add_patch(FancyArrowPatch((0.4,0.5),(1.25,0.5),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(0.8,0.5,'Airflow\nDirection',ha='center',va='center')

pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.savefig('images/HXTubesTerms.pdf')
pylab.show()