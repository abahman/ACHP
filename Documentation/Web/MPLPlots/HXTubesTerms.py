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

#plot lower tube
verts1 = [
    (0.24, 0.05),  # P0: draw the outer tube line
    (0.3, 0.05), # P1
    (0.3, -0.05), # P2
    (0.24, -0.05), # P3
    (-0.24, -0.05), # P4
    (-0.3, -0.05), # P5
    (-0.3, 0.05), # P6
    (-0.24, 0.05), ## P7
    (0.24, 0.05),  # end: close the loop
    (0.22, 0.035),  # P0 : start of the first port on the right
    (0.28, 0.035), # P1
    (0.28, -0.035), # P2
    (0.22, -0.035), # P3
    (0.15, -0.035), # P4
    (0.15, 0.035), # P5
    (0.22, 0.035), # closepoly 
    (-0.22, -0.035), # P0 : draw the port on the left
    (-0.28, -0.035), # P1
    (-0.28, 0.035), # P2
    (-0.22, 0.035), # P3
    (-0.15,0.035), # P4
    (-0.15,-0.035), # P5
    (-0.22,-0.035),  # closepoly
    (0.13, 0.035), # draw the rectangle on the right
    (0.13, -0.035),
    (0.01,-0.035),
    (0.01,0.035),
    (0.13,0.035),
    (-0.13, 0.035), #draw the revtangle on the left
    (-0.13, -0.035),
    (-0.01,-0.035),
    (-0.01,0.035),
    (-0.13,0.035),
    ]
#code of plot lower tube
codes1 = [Path.MOVETO, # draw the outter line tube
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY, #end
         Path.MOVETO, #first port on the right
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end of plot
         Path.MOVETO, #start ro draw the port on the left
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end of plot
         Path.MOVETO, #start to draw the rectangle on the right
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end
         Path.MOVETO, # start to plot rectangle on the left
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end
         ]

#plot upper tube
verts2 = [
    (0.24, 0.555),  # P0 : draw the outer tube of the upper tube
    (0.3, 0.555), # P1
    (0.3, 0.455), # P2
    (0.24, 0.455), # P3
    (-0.24, 0.455), # P4
    (-0.3, 0.455), # P5
    (-0.3, 0.555), # P6
    (-0.24, 0.555), ## P7
    (0.24, 0.555),  # to close the loop , please ignored because it is the same as the first point
    (0.22, 0.54),  # P0 : start of the first port on the right
    (0.28, 0.54), # P1
    (0.28, 0.47), # P2
    (0.22, 0.47), # P3
    (0.15, 0.47), # P4
    (0.15, 0.54), # P5
    (0.22, 0.54), # closepoly 
    (-0.22, 0.47), # P0 : draw the port on the left
    (-0.28, 0.47), # P1
    (-0.28, 0.54), # P2
    (-0.22, 0.54), # P3
    (-0.15,0.54), # P4
    (-0.15,0.47), # P5
    (-0.22,0.47),  # closepoly
    (0.13, 0.54), # draw the rectangle on the right
    (0.13, 0.47),
    (0.01, 0.47),
    (0.01, 0.54),
    (0.13, 0.54),
    (-0.13, 0.54), #draw the revtangle on the left
    (-0.13, 0.47),
    (-0.01, 0.47),
    (-0.01, 0.54),
    (-0.13, 0.54), 
    ]

#code of plot upper tube
codes2 = [Path.MOVETO, #draw the outer tube of the upper tube
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         Path.MOVETO, #first port on the right
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end of plot
         Path.MOVETO, #start to draw the port on the left
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end of plot
         Path.MOVETO, #start to draw the rectangle on the right
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end
         Path.MOVETO, # start to plot rectangle on the left
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY, #end
         ]

#plot the red fin
verts3 = [
    (0.28, 0.45),  # P0 : draw the red rectangle
    (0.28, 0.055), # P1
    (-0.28, 0.055), # P2
    (-0.28, 0.45), # P3
    (0.28, 0.45), # close loop
    
    (0.2, 0.4),  # P0 : draw the 1st vertical louver fins from right
    (0.2, 0.41), # P1
    (0.22, 0.41), # P2
    (0.22, 0.4), # P3
    (0.22, 0.1), # P4
    (0.22, 0.09), # P5
    (0.2, 0.09), # P6
    (0.2, 0.1), ## P7
    (0.2, 0.4), #close loop
    
    (0.16, 0.4),  # P0 : draw the 2nd vertical louver fins from right
    (0.16, 0.41), # P1
    (0.14, 0.41), # P2
    (0.14, 0.4), # P3
    (0.14, 0.1), # P4
    (0.14, 0.09), # P5
    (0.16, 0.09), # P6
    (0.16, 0.1), ## P7
    (0.16, 0.4), #close loop
    
    (0.1, 0.4),  # P0 : draw the 3rd vertical louver fins from right
    (0.1, 0.41), # P1
    (0.08, 0.41), # P2
    (0.08, 0.4), # P3
    (0.08, 0.1), # P4
    (0.08, 0.09), # P5
    (0.1, 0.09), # P6
    (0.1, 0.1), ## P7
    (0.1, 0.4), #close loop
    
    (-0.2, 0.4),  # P0 : draw the 1st vertical louver fins from left
    (-0.2, 0.41), # P1
    (-0.22, 0.41), # P2
    (-0.22, 0.4), # P3
    (-0.22, 0.1), # P4
    (-0.22, 0.09), # P5
    (-0.2, 0.09), # P6
    (-0.2, 0.1), ## P7
    (-0.2, 0.4), #close loop
    
    (-0.16, 0.4),  # P0 : draw the 2nd vertical louver fins from left
    (-0.16, 0.41), # P1
    (-0.14, 0.41), # P2
    (-0.14, 0.4), # P3
    (-0.14, 0.1), # P4
    (-0.14, 0.09), # P5
    (-0.16, 0.09), # P6
    (-0.16, 0.1), ## P7
    (-0.16, 0.4), #close loop
    
    (-0.1, 0.4),  # P0 : draw the 3rd vertical louver fins from left
    (-0.1, 0.41), # P1
    (-0.08, 0.41), # P2
    (-0.08, 0.4), # P3
    (-0.08, 0.1), # P4
    (-0.08, 0.09), # P5
    (-0.1, 0.09), # P6
    (-0.1, 0.1), ## P7
    (-0.1, 0.4), #close loop
    ]

#code to plot red fin
codes3 = [Path.MOVETO, #draw the outer tube of the upper tube
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         
         Path.MOVETO, #draw first vertical louver fins from right
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         
         Path.MOVETO, #draw 2nd vertical louver fins from right
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         
         Path.MOVETO, #draw 3rd vertical louver fins from right
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         
         Path.MOVETO, #draw first vertical louver fins from left
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         
         Path.MOVETO, #draw 2nd vertical louver fins from left
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         
         Path.MOVETO, #draw 3rd vertical louver fins from left
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.LINETO,
         Path.CURVE4,
         Path.CURVE4,
         Path.CURVE4,
         Path.CLOSEPOLY,
         ]

path1 = Path(verts1, codes1)
path2 = Path(verts2, codes2)
path3 = Path(verts3, codes3)
ax = fig.add_subplot(111)
patch1 = patches.PathPatch(path1, facecolor='none',edgecolor='blue', lw=0.8)
patch2 = patches.PathPatch(path2, facecolor='none',edgecolor='blue', lw=0.8)
patch3 = patches.PathPatch(path3, facecolor='none',edgecolor='red', lw=0.8)

ax.add_patch(patch1)
ax.add_patch(patch2)
ax.add_patch(patch3)

#Uncomment the next two line to see the point's path used
#xs, ys = zip(*verts1)
#pylab.plot(xs, ys, 'x--', lw=0.8, color='black', ms=10)

##Dimension lines
pylab.plot(np.r_[-.3,-.3],np.r_[0.5,0.95],'k')
pylab.plot(np.r_[0.3,0.3],np.r_[0.5,0.95],'k')
pylab.text(0,0.85,'Tube\nDepth\n$T_d$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((-0.3,0.7),(0.3,0.7),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.plot(np.r_[-0.7,-0.3],np.r_[0.45,0.45],'k')
pylab.plot(np.r_[-0.7,-0.3],np.r_[0.05,0.05],'k')
pylab.text(-0.75,0.25,'Tube\nSpacing\n$b$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((-0.5,0.05),(-0.5,0.45),arrowstyle='<|-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

pylab.plot(np.r_[0.5,0.3],np.r_[0.05,0.05],'k')
pylab.plot(np.r_[0.5,0.3],np.r_[-0.05,-0.05],'k')
pylab.text(0.7,0.0,'Major\nDiameter\n$H_t$',ha='center',va='center')
pylab.gca().add_patch(FancyArrowPatch((0.4,0.05),(0.4,0.2),arrowstyle='<|-',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.gca().add_patch(FancyArrowPatch((0.4,-0.2),(0.4,-0.05),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))

#Airflow arrow
pylab.gca().add_patch(FancyArrowPatch((0.4,0.25),(1.25,0.25),arrowstyle='-|>',fc='k',ec='k',mutation_scale=20,lw=0.8))
pylab.text(0.8,0.25,'Airflow\nDirection',ha='center',va='center')

pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.savefig('images/HXTubesTerms.pdf')
pylab.show()