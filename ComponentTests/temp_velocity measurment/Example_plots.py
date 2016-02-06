'''
Created on Jan 4, 2016

@author: ammarbahman

Note: this file plots the countor of velocity profile in 60K ECU

'''

#===============================================================================
# EXAMPLES
#===============================================================================

# from pylab import meshgrid, sin, cos, linspace, contourf, savefig, clf
# import matplotlib.pyplot as plt
# #plt.rc('font',family='serif',size =11.0)
# plt.rc('text', usetex=True)
# 
# x, y = meshgrid(*(linspace(-1,1,500),)*2)
# z = sin(20*x**2)*cos(30*y)
# c = contourf(x,y,z,50)
# savefig('full_vector.pdf')
# clf()
# c = contourf(x,y,z,50,rasterized=True)
# savefig('rasterized.pdf')


# import pylab
# from pylab import arange,pi,sin,cos,sqrt
# fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
# inches_per_pt = 1.0/72.27               # Convert pt to inch
# golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
# fig_width = fig_width_pt*inches_per_pt  # width in inches
# fig_height = fig_width*golden_mean      # height in inches
# fig_size =  [fig_width,fig_height]
# params = {'backend': 'ps',
#           'axes.labelsize': 10,
#           'text.fontsize': 10,
#           'legend.fontsize': 10,
#           'xtick.labelsize': 8,
#           'ytick.labelsize': 8,
#           'text.usetex': True,
#           'figure.figsize': fig_size}
# pylab.rcParams.update(params)
# # Generate data
# x = pylab.arange(-2*pi,2*pi,0.01)
# y1 = sin(x)
# y2 = cos(x)
# # Plot data
# pylab.figure(1)
# pylab.clf()
# pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
# pylab.plot(x,y1,'g:',label='$\sin(x)$')
# pylab.plot(x,y2,'-b',label='$\cos(x)$')
# pylab.xlabel('$x$ (radians)')
# pylab.ylabel('$y$')
# pylab.legend()
# pylab.savefig('fig1.pdf')


# import pylab
# from pylab import arange,pi,sin,cos,sqrt
# fig_width_pt = 246.0  # Get this from LaTeX using \showthe\columnwidth
# inches_per_pt = 1.0/72.27               # Convert pt to inch
# golden_mean = (sqrt(5)-1.0)/2.0         # Aesthetic ratio
# fig_width = fig_width_pt*inches_per_pt  # width in inches
# fig_height = fig_width*golden_mean      # height in inches
# fig_size =  [fig_width,fig_height]
# params = {'backend': 'ps',
#           'axes.labelsize': 10,
#           'text.fontsize': 10,
#           'legend.fontsize': 10,
#           'xtick.labelsize': 8,
#           'ytick.labelsize': 8,
#           'text.usetex': True,
#           'figure.figsize': fig_size}
# pylab.rcParams.update(params)
# # Generate data
# x = pylab.arange(-2*pi,2*pi,0.01)
# y1 = sin(x)
# y2 = cos(x)
# # Plot data
# # Plot data
# pylab.figure(1)
# pylab.clf()
# pylab.axes([0.125,0.2,0.95-0.125,0.95-0.2])
# pylab.plot(x,y1,'g:',label=r'sin(x)')
# pylab.plot(x,y2,'-b',label=r'cos(x)')
# pylab.xlabel(r'x (radians)')
# pylab.ylabel(r'y')
# pylab.legend()
# pylab.savefig('fig2.pdf')


# from matplotlib import rc
# import numpy as np
# from pylab import figure, axes, plot, xlabel, ylabel, title, grid, savefig, show
# 
# rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Roman']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)
# figure(1)
# ax = axes([0.1, 0.1, 0.8, 0.7])
# t = np.arange(0.0, 1.0+0.01, 0.01)
# s = np.cos(2*2*np.pi*t)+2
# plot(t, s)
# 
# xlabel(r'\textbf{time (s)}')
# ylabel(r'\textit{voltage (mV)}',fontsize=16)
# title(r"\TeX\ is Number $\displaystyle\sum_{n=1}^\infty\frac{-e^{i\pi}}{2^n}$!", fontsize=16, color='r')
# grid(True)
# savefig('tex_demo.pdf')
# show()


# import numpy as np
# import matplotlib.cm as cm
# import matplotlib.mlab as mlab
# import matplotlib.pyplot as plt
# 
# delta = 0.025
# x = y = np.arange(-3.0, 3.0, delta)
# X, Y = np.meshgrid(x, y)
# Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
# Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
# Z = Z2 - Z1  # difference of Gaussians
# 
# im = plt.imshow(Z, interpolation='bilinear', cmap=cm.RdYlGn,
#                 origin='lower', extent=[-3, 3, -3, 3],
#                 vmax=abs(Z).max(), vmin=-abs(Z).max())
# plt.savefig('countor.pdf')
# plt.show()


# from scipy.interpolate import griddata
# import matplotlib.pyplot as plt
# np.random.seed(0)
# 
# x = np.random.normal(size=200)
# y = np.random.normal(size=200)
# v = np.sqrt(x**2+y**2)
# 
# xg, yg = np.mgrid[x.min():x.max():100j, y.min():y.max():100j]
# vg = griddata((x, y), v, (xg, yg), method='cubic')
# plt.contourf(xg, yg, vg)
# plt.scatter(x, y, c=v)
# plt.savefig('countor_scatter.pdf')
# plt.show()

# import pylab as plt
# 
# x=[12, 13, 14, 15, 16] # x-axis coordinates
# y=[14, 15, 16, 17, 18] # y-axis coordinates
# 
# v_x=[6, 6, 6, 6, 6] # x-component of velocity
# v_y=[1,4,3,2,1]     # y-component of velocity
# 
# plt.quiver(x,y,v_x,v_y)
# plt.xlim(11,17)
# plt.ylim(13,19)
# 
# plt.savefig('quiver.pdf')
# plt.show()

#===============================================================================
# END OF EXMAPLES
#===============================================================================