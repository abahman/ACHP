'''
Created on Jan 4, 2016

@author: ammarbahman

Note: this file plots the countor of velocity profile in 60K ECU

'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from pylab import contourf, clf

import matplotlib as mpl
#mpl.use('pgf')

def figsize(scale):
    fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
    inches_per_pt = 1.0/72.27                       # Convert pt to inch
    golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
    fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
    fig_height = fig_width*golden_mean              # height in inches
    fig_size = [fig_width,fig_height]
    return fig_size

pgf_with_latex = {                      # setup matplotlib to use latex for output
"pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
"text.usetex": True,                # use LaTeX to write all text
"font.family": "serif",
"font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
"font.sans-serif": [],
"font.monospace": [],
"axes.labelsize": 10,               # LaTeX default is 10pt font.
"font.size": 10,
"legend.fontsize": 8,               # Make the legend/label fonts a little smaller
"xtick.labelsize": 8,
"ytick.labelsize": 8,
"figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
"pgf.preamble": [
r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
        ]
    }
mpl.rcParams.update(pgf_with_latex)

 
x = np.array([6.219,12.438,18.657])
y = np.array([20.625,16.875,13.125,9.375,5.625,1.875])

X, Y = np.meshgrid(x, y)
#tempearture profile for Test 5
T = np.matrix([[12.45,12.45,12.91],[13.69,15.29,14.65],[10.99,10.97,10.59],[10.76,10.72,
             10.58],[7.849,8.311,8.868],[4.798,5.016,6.926]])

im = plt.imshow(T, interpolation='bicubic',extent=[0, 24.875, 0, 22.5],
                vmax=18, vmin=0)

cbar = plt.colorbar(im)
cbar.ax.set_ylabel(r'Temperature [\textdegree C]')
plt.ylim(0,22.5)
plt.xlim(0,24.875)
plt.xticks([0, 5, 10, 15, 20, 24.875],
          [r'$0$', r'$5$', r'$10$', r'$15$', r'$20$', r'$24.875$'])
plt.yticks([0, 5, 10, 15, 20, 22.5],
          [r'$0$', r'$5$', r'$10$', r'$15$', r'$20$', r'$22.5$'])
plt.xlabel('Evaporator width [in]')
plt.ylabel('Evaporator height [in]')
#plt.title('Temperature profile of Test 5')


### TO SHOW the values with the measurment grid on the plot
# for i in range(len(x)):
#     for j in range(len(y)):
#         plt.plot(x[i],y[j],'ko')
#         plt.annotate(T[j,i], (x[i],y[j]))
    
plt.savefig('temp_profile/temp_profile_test5.pdf')
plt.show()



plt.figure()
#tempearture profile for Test 6
T = np.matrix([[12.14,11.36,11.1],[13.17,13.61,14.02],[10.4,11.43,11.51],
               [9.983,10.64,11.01],[7.872,7.648,8.034],[5.318,5.505,5.93]])

im = plt.imshow(T, interpolation='bicubic',extent=[0, 24.875, 0, 22.5],
                vmax=18, vmin=0)

cbar = plt.colorbar(im)
cbar.ax.set_ylabel(r'Temperature [\textdegree C]')
plt.ylim(0,22.5)
plt.xlim(0,24.875)
plt.xticks([0, 5, 10, 15, 20, 24.875],
          [r'$0$', r'$5$', r'$10$', r'$15$', r'$20$', r'$24.875$'])
plt.yticks([0, 5, 10, 15, 20, 22.5],
          [r'$0$', r'$5$', r'$10$', r'$15$', r'$20$', r'$22.5$'])
plt.xlabel('Evaporator width [in]')
plt.ylabel('Evaporator height [in]')
#plt.title('Temperature profile of Test 6')


### TO SHOW the values with the measurment grid on the plot
# for i in range(len(x)):
#     for j in range(len(y)):
#         plt.plot(x[i],y[j],'ko')
#         plt.annotate(T[j,i], (x[i],y[j]))
    
plt.savefig('temp_profile/temp_profile_test6.pdf')
plt.show()