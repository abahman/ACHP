import pylab,numpy as np
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import matplotlib.patches as patches
import matplotlib.pyplot as plt

fig=pylab.figure(figsize=(6,5))


#plot the fin
verts1 = [
    (0.3, 0.45),  # P0 : draw the red rectangle
    (0.3, 0.055), # P1
    (-0.3, 0.055), # P2
    (-0.3, 0.45), # P3
    (0.3, 0.45), # close loop
    
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

#code to plot regtangle fin
codes1 = [
         Path.MOVETO, #draw the outer tube of the upper tube
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


verts2 = [
    (0.3, -0.5),  # draw the cross section louver fin
    (0.3, -0.49),  
    (0.25, -0.49), 
    (0.23, -0.45),
    (0.22, -0.45), 
    (0.245, -0.5),
    (0.3, -0.5),
    
    (0.22, -0.5), #move to : start 2nd louver
    (0.195, -0.45),
    (0.185,-0.45),
    (0.21,-0.5),
    (0.22, -0.5),
   
    (0.18, -0.5), #move to : start 3nd louver
    (0.155, -0.45),
    (0.145,-0.45),
    (0.17,-0.5),
    (0.18, -0.5),
    
    (0.14, -0.5), #move to : start center louver
    (0.115, -0.46),
    (0.0,-0.46),
    (0.14,-0.5), #move to start point
    (0.13, -0.5),
    (0.11,-0.47),
    (0.0,-0.47),
    
    #REPEAT FOR THE REFLECTION ON THE LEFT HAND SIDE
    (-0.3, -0.5),  # draw the cross section louver fin
    (-0.3, -0.49),  
    (-0.25, -0.49), 
    (-0.23, -0.45),
    (-0.22, -0.45), 
    (-0.245, -0.5),
    (-0.3, -0.5),
    
    (-0.22, -0.5), #move to : start 2nd louver
    (-0.195, -0.45),
    (-0.185,-0.45),
    (-0.21,-0.5),
    (-0.22, -0.5),
   
    (-0.18, -0.5), #move to : start 3nd louver
    (-0.155, -0.45),
    (-0.145,-0.45),
    (-0.17,-0.5),
    (-0.18, -0.5),
    
    (-0.14, -0.5), #move to : start center louver
    (-0.115, -0.46),
    (0.0,-0.46),
    (-0.14,-0.5), #move to start point
    (-0.13, -0.5),
    (-0.11,-0.47),
    (0.0,-0.47),
    
    ]

#code to plot louvered fin
codes2 = [
         Path.MOVETO, # draw the cross section fin
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         
         Path.MOVETO, ##move to : start 2nd louver
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,

         Path.MOVETO, ##move to : start 3rd louver
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         
         Path.MOVETO, ##move to : start 3rd louver
         Path.LINETO,
         Path.LINETO,
         Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         
         #REPEAT FOR THE REFLECTION ON THE LEFT HAND SIDE
         Path.MOVETO, # draw the cross section fin
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         
         Path.MOVETO, ##move to : start 2nd louver
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,

         Path.MOVETO, ##move to : start 3rd louver
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         Path.CLOSEPOLY,
         
         Path.MOVETO, ##move to : start 3rd louver
         Path.LINETO,
         Path.LINETO,
         Path.MOVETO,
         Path.LINETO,
         Path.LINETO,
         Path.LINETO,
         ]

path1 = Path(verts1, codes1)
path2 = Path(verts2, codes2)
ax = fig.add_subplot(111)
patch1 = patches.PathPatch(path1, facecolor='none',edgecolor='black', lw=0.8)
patch2 = patches.PathPatch(path2, facecolor='none',edgecolor='black', lw=0.8)

ax.add_patch(patch1)
ax.add_patch(patch2)



#h=0.45
#x=np.array([0,1,2,3,4])
#y=np.array([0,h,0,h,0])
# pylab.plot(x,y,lw=12,color='grey')
# pylab.plot(x,y+2*h,lw=12,color='grey')
# pylab.plot(x,y+4*h,lw=12,color='grey')

#centerlines of each fin
# pylab.plot(x,y,'--',lw=1,color='k')
# pylab.plot(x,y+2*h,'--',lw=1,color='k')
# pylab.plot(x,y+4*h,'--',lw=1,color='k')

#Label for Llouv
h=0.41
xx=-0.36
pylab.plot(np.r_[xx,xx],np.r_[0.09,h],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[0.09,0.09],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[h,h],'b')
pylab.text(xx-0.01,0.5*(h+0.09),'$L_{louv}$',ha='right',va='center')

#Label for sf
h=0.45
xx=0.36
pylab.plot(np.r_[xx,xx],np.r_[0.055,h],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[0.055,0.055],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[h,h],'b')
pylab.text(xx+0.01,0.5*(h),'$s_f$',ha='left',va='center')

# #Label for pd
# xx=-0.25
# pylab.plot(np.r_[xx,xx],np.r_[4*h,5*h],'b')
# pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[4*h,4*h],'b')
# pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[5*h,5*h],'b')
# pylab.text(xx,4.5*h,'$p_d$',ha='right',va='center')
# 
#Label for Lf
yy=0.51
pylab.plot(np.r_[-0.3,0.3],np.r_[yy,yy],'b')
pylab.plot(np.r_[-0.3,-0.3],np.r_[yy-0.05,yy+0.05],'b')
pylab.plot(np.r_[0.3,0.3],np.r_[yy-0.05,yy+0.05],'b')
pylab.text(0.0,yy+0.01,'$L_f$',ha='center',va='bottom')
#
#Label for lp
yy=-0.01
pylab.plot(np.r_[0.2,0.22],np.r_[yy,yy],'b')
pylab.plot(np.r_[0.2,0.2],np.r_[yy-0.05,yy+0.05],'b')
pylab.plot(np.r_[0.22,0.22],np.r_[yy-0.05,yy+0.05],'b')
pylab.text(0.17,yy-0.04,'$l_p$',ha='center',va='bottom')
#
#Label for Lalpha
xxx=-0.22
h=0.41
yy=-0.5
#pylab.plot(np.r_[xxx,xxx],np.r_[0.1,yy],'b')
pylab.plot(np.r_[xxx-0.02,xxx+0.016],np.r_[yy-0.01,yy+0.06],'b')
pylab.plot(np.r_[xxx-0.02,xxx+0.05],np.r_[yy-0.01,yy-0.01],'b')
pylab.text(xxx+0.05,yy-0.07,'$L_{alpha}$',ha='right',va='center')
#
#Label for delta
h=-0.5
xx=0.36
pylab.plot(np.r_[xx,xx],np.r_[h,h+0.01],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[h,h],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[h+0.01,h+0.01],'b')
pylab.text(xx+0.01,h+0.04,'$\delta$',ha='left',va='center') 
#
#Label for lh
h=-0.5
xx=-0.36
pylab.plot(np.r_[xx,xx],np.r_[h,h+0.05],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[h,h],'b')
pylab.plot(np.r_[xx-0.05,xx+0.05],np.r_[h+0.05,h+0.05],'b')
pylab.text(xx-0.05,h+0.025,'$l_h$',ha='right',va='center') 
# #Label for s
# xxx=4
# pylab.plot(np.r_[xxx,xxx],np.r_[0.1,2*(h-0.1)],'b')
# pylab.plot(np.r_[xxx-0.05,xxx+0.05],np.r_[0.1,0.1],'b')
# pylab.plot(np.r_[xxx-0.05,xxx+0.05],np.r_[2*(h-0.1),2*(h-0.1)],'b')
# pylab.text(xxx,h,'$s$',ha='right',va='center')
# 
# 
# #Label for t - would probably be better to do this with arrows, but so be it
# pylab.plot(np.r_[2,2.6],np.r_[4*h-0.08,4*h-0.08],'b')
# pylab.plot(np.r_[2,2.6],np.r_[4*h+0.1,4*h+0.1],'b')
# pylab.text(2.6,4*h+0.01,'$t$',va='center',ha='left')

pylab.gca().axis('equal')
pylab.gca().axis('off')
pylab.savefig('images/FinTermDefinition.pdf')
pylab.show()