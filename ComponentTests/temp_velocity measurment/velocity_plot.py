'''
Created on Jan 4, 2016

@author: ammarbahman

Note: this file plots the countor of temperature profile in 60K ECU

'''
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import matplotlib.mlab as mlab
import matplotlib.cm as cm
from pylab import contourf, clf

#===============================================================================
# Latex render
#===============================================================================
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
#===============================================================================
# END of Latex render
#===============================================================================
 


x = np.array([1.517,5.846,10.870,14.005,19.029,23.358])
y = np.array([21.308,17.933,14.265,11.250,8.235,4.568,1.193])

X, Y = np.meshgrid(x, y)
#===============================================================================
# Velocity profile for All tests
#===============================================================================
Test =['4','5','Fan only']

V4_max = np.matrix([[3.10,2.49,5.62,3.61,4.43,6.52],
                    [2.83,3.43,3.61,2.48,4.62,3.28],
                    [3.74,4.28,4.38,3.37,3.82,3.47],
                    [4.41,4.18,4.56,4.52,4.29,4.70],
                    [4.01,4.81,5.11,4.06,4.51,4.32],
                    [4.24,4.12,3.90,4.40,4.25,4.25],
                    [3.22,4.21,4.73,4.38,5.68,3.84]])
V4_min = np.matrix([[2.56,2.31,5.15,2.46,4.17,5.43],
                    [1.92,2.85,3.27,2.20,2.27,1.86],
                    [1.40,3.80,3.79,3.11,2.99,2.40],
                    [2.08,3.53,3.74,3.80,2.98,3.40],
                    [3.27,3.85,4.16,3.69,4.12,3.49],
                    [3.66,3.15,3.08,3.43,3.96,3.45],
                    [2.12,3.82,3.56,4.09,5.25,3.44]])

V5_max = np.matrix([[1.77,5.05,3.49,1.12,3.66,5.26],
                    [3.11,2.55,2.61,2.39,2.88,3.65],
                    [2.48,2.49,3.58,1.91,2.58,1.90],
                    [4.35,5.09,4.94,4.84,4.80,4.82],
                    [3.43,3.43,3.63,4.42,3.65,2.56],
                    [3.54,3.08,3.98,3.87,2.81,2.29],
                    [3.84,3.75,3.67,4.14,4.42,3.46]])
V5_min = np.matrix([[1.58,4.78,3.13,1.01,3.30,4.76],
                    [2.78,2.38,2.49,2.31,2.42,3.44],
                    [2.30,2.26,3.43,1.70,2.41,1.79],
                    [4.25,4.99,4.84,4.72,4.69,4.59],
                    [3.13,3.21,3.41,4.09,3.34,2.36],
                    [3.32,2.88,3.66,3.69,2.66,2.20],
                    [3.33,3.21,3.45,3.93,4.19,3.23]])

VFan_max = np.matrix([[3.18,4.32,2.14,3.11,3.93,4.51],
                      [3.60,2.36,1.77,2.12,1.56,2.14],
                      [1.76,2.94,2.28,3.36,2.71,2.06],
                      [4.04,4.78,4.97,4.72,4.75,4.28],
                      [2.89,3.69,4.15,3.78,4.14,3.18],
                      [1.86,2.86,3.11,3.16,3.37,2.56],
                      [2.86,2.70,3.24,3.38,3.24,2.55]])
VFan_min = np.matrix([[2.90,4.15,2.04,2.73,3.78,4.44],
                      [3.47,2.09,1.61,1.97,1.41,2.00],
                      [1.63,2.77,2.20,3.26,2.62,1.91],
                      [3.92,4.68,4.86,4.62,4.63,4.18],
                      [2.74,3.33,4.04,3.65,4.03,3.05],
                      [1.75,2.66,2.92,3.05,3.25,2.46],
                      [2.65,2.62,3.12,3.31,3.15,2.48]])

#Average the data
V4 = (V4_max+V4_min)/2
V5 = (V5_max+V5_min)/2
VFan = (VFan_max+VFan_min)/2

V_data = [V4,V5,VFan]

for i in range(len(Test)):
    plt.figure()
    im = plt.imshow(V_data[i], interpolation='bicubic',extent=[0, 24.875, 0, 22.5],
                vmax=6, vmin=0)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Velocity [m/s]')
    plt.ylim(0,22.5)
    plt.xlim(0,24.875)
    plt.xticks([0, 5, 10, 15, 20, 24.875],
              [r'$0$', r'$5$', r'$10$', r'$15$', r'$20$', r'$24.875$'])
    plt.yticks([0, 4, 8, 12, 16, 20, 22.5],
              [r'$0$', r'$4$', r'$8$', r'$12$',r'$16$', r'$20$', r'$22.5$'])
    plt.xlabel('Evaporator width [in]')
    plt.ylabel('Evaporator height [in]')
    plt.title('Velocity profile of Test '+Test[i])
      
 
    ### TO SHOW the values with the measurment grid on the plot
#     for k in range(len(x)):
#         for j in range(len(y)):
#             plt.plot(x[k],y[j],'ko')
#             plt.annotate(V_data[i][j,k], (x[k],y[j]))
          
    plt.savefig('velocity_profile/velocity_profile_test'+Test[i]+'.pdf')
    plt.show()

#===============================================================================
# TO SHOW all plots in one Figure
#===============================================================================
fig = plt.figure(1, figsize=(12.5, 3))
for i in range(len(Test)):
    ax = plt.subplot(1, 3, i+1)
    im = plt.imshow(V_data[i], interpolation='bicubic',extent=[0, 24.875, 0, 22.5],
                vmax=6, vmin=0)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel(r'Velocity [m/s]')
    plt.ylim(0,22.5)
    plt.xlim(0,24.875)
    plt.xticks([0, 5, 10, 15, 20, 24.875],
              [r'$0$', r'$5$', r'$10$', r'$15$', r'$20$', r'$24.875$'])
    plt.yticks([0, 4, 8, 12, 16, 20, 22.5],
              [r'$0$', r'$4$', r'$8$', r'$12$',r'$16$', r'$20$', r'$22.5$'])
    plt.xlabel('Evaporator width [in]')
    plt.ylabel('Evaporator height [in]')
    plt.title('Velocity profile of Test '+Test[i])
   
        
    ### TO SHOW the values with the measurment grid on the plot
    # for k in range(len(x)):
    #     for j in range(len(y)):
    #         plt.plot(x[k],y[j],'ko')
    #         plt.annotate(V_data[i][j,k], (x[k],y[j]))
fig.set_tight_layout(True)
plt.savefig('velocity_profile/velocity_profile_combined.pdf')
plt.show()