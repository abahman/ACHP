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


import numpy as np
import matplotlib.cm as cm
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
 
delta = 0.025
x = y = np.arange(-3.0, 3.0, delta)
X, Y = np.meshgrid(x, y)
Z1 = mlab.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = mlab.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
Z = Z2 - Z1  # difference of Gaussians
 
im = plt.imshow(Z, interpolation='bilinear', cmap=cm.RdYlGn,
                origin='lower', extent=[-3, 3, -3, 3],
                vmax=abs(Z).max(), vmin=-abs(Z).max())
plt.savefig('countor.pdf')
plt.show()
 
 
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
np.random.seed(0)
 
x = np.random.normal(size=200)
y = np.random.normal(size=200)
v = np.sqrt(x**2+y**2)
 
xg, yg = np.mgrid[x.min():x.max():100j, y.min():y.max():100j]
vg = griddata((x, y), v, (xg, yg), method='cubic')
plt.contourf(xg, yg, vg)
plt.scatter(x, y, c=v)
plt.savefig('countor_scatter.pdf')
plt.show()

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

# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# from matplotlib.ticker import LinearLocator, FormatStrFormatter
# import matplotlib.pyplot as plt
# import numpy as np
# 
# fig = plt.figure()
# ax = fig.gca(projection='3d')
# x = np.array([6.219,12.438,18.657])
# y = np.array([20.625,16.875,13.125,9.375,5.625,1.875])
# X, Y = np.meshgrid(x, y)
# Z = np.matrix([[15.68,15.79,14.86],[15.64,16.08,15.95],[14.6,14.83,14.57],
#                 [14.7,14.86,14.65],[13.63,13.67,13.75],[11.92,12.3,12.6]])
# 
# surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
#         linewidth=0, antialiased=False)
# ax.set_zlim(4, 22)
# 
# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
# 
# fig.colorbar(surf)
# 
# plt.show()


import scipy
import matplotlib.pyplot as plt
import numpy as np

x = scipy.array([ 21.308,17.933,14.265,11.25,8.235,4.568,1.193])
y = scipy.array([3.9875,2.885,3.37916667,3.84916667,4.11666667,3.82416667,4.02833333])
result = scipy.poly1d([0.0]) #setting result = 0

for i in range(0,len(x)): #number of polynomials L_k(x).
    temp_numerator = scipy.poly1d([1.0]) # resets temp_numerator such that a new numerator can be created for each i.
    denumerator = 1.0 #resets denumerator such that a new denumerator can be created for each i.
    for j in range(0,len(x)):
        if i != j:
            temp_numerator *= scipy.poly1d([1.0,-x[j]]) #finds numerator for L_i
            denumerator *= x[i]-x[j] #finds denumerator for L_i
    result += (temp_numerator/denumerator) * y[i] #linear combination

print("The result is: ")
print(result)

x_val = np.arange(min(x),max(x)+1, 0.1) #generates x values we would like to evaluate.
plt.xlabel('x'); plt.ylabel('p(x)')
plt.grid(True)
for i in range(0,len(x)):
    plt.plot([x[i]],[y[i]],'ro') #plot the points
plt.plot(x_val, result(x_val)) #result(x_val) gives the value of our Lagrange polynomial.
plt.axis([min(x)-1, max(x)+1, min(y)-1, max(y)+1])
plt.show()

#===============================================================================
# END OF EXMAPLES
#===============================================================================