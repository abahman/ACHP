'''
Created on March 7, 2017

@author: ammarbahman

Note: this file plots the superheat  in 60K ECU for IJR paper

'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata, interp1d
from scipy import exp
import scipy
import scipy.ndimage as ndimage
#import matplotlib.mlab as mlab
#import matplotlib.cm as cm
#from pylab import contourf, clf
import matplotlib as mpl
mpl.style.use('classic')

#===============================================================================
# Latex render
#===============================================================================
import matplotlib as mpl
from numpy import integer
from numba.targets.randomimpl import f_impl
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
"legend.labelspacing":0.2,
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
 

#import the excel file
df = pd.read_excel("superheat/superheat.xlsx")
B_test1 = np.array(df[:6]["Test1"])
M_test1 = np.array(df[10:16]["Test1"])
I_test1 = np.array(df[20:26]["Test1"])

B_test2 = np.array(df[:6]["Test2"])
M_test2 = np.array(df[10:16]["Test2"])
I_test2 = np.array(df[20:26]["Test2"])

B_test3 = np.array(df[:6]["Test3"])
M_test3 = np.array(df[10:16]["Test3"])
I_test3 = np.array(df[20:26]["Test3"])

B_test4 = np.array(df[:6]["Test4"])
M_test4 = np.array(df[10:16]["Test4"])
I_test4 = np.array(df[20:26]["Test4"])

B_test5 = np.array(df[:6]["Test5"])
M_test5 = np.array(df[10:16]["Test5"])
I_test5 = np.array(df[20:26]["Test5"])

B_test6 = np.array(df[:6]["Test6"])
M_test6 = np.array(df[10:16]["Test6"])
I_test6 = np.array(df[20:26]["Test6"])

B_testB = np.array(df[:6]["TestB"])
M_testB = np.array(df[10:16]["TestB"])
I_testB = np.array(df[20:26]["TestB"])

B_testC = np.array(df[:6]["TestC"])
M_testC = np.array(df[10:16]["TestC"])
I_testC = np.array(df[20:26]["TestC"])

    
    
##########plot superheat -- Baseline##########
plt.plot(np.arange(1,7,1),B_testC,'--ko',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test C')
plt.errorbar(np.arange(1,7,1),B_testC,yerr=0.5,fmt='',linestyle="None",color='k')
plt.plot(np.arange(1,7,1),B_testB,'--bs',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test B')
plt.errorbar(np.arange(1,7,1),B_testB,yerr=0.5,fmt='',linestyle="None",color='b')
plt.plot(np.arange(1,7,1),B_test6,'--r^',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 6')
plt.errorbar(np.arange(1,7,1),B_test6,yerr=0.5,fmt='',linestyle="None",color='r')
plt.plot(np.arange(1,7,1),B_test5,'--g*',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 5')
plt.errorbar(np.arange(1,7,1),B_test5,yerr=0.5,fmt='',linestyle="None",color='g')
plt.plot(np.arange(1,7,1),B_test4,'--P',markersize=5,markeredgewidth=0.1,alpha=0.9,color='brown',label=r'Test 4')
plt.errorbar(np.arange(1,7,1),B_test4,yerr=0.5,fmt='',linestyle="None",color='brown')
plt.plot(np.arange(1,7,1),B_test3,'--cH',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 3')
plt.errorbar(np.arange(1,7,1),B_test3,yerr=0.5,fmt='',linestyle="None",color='c')
plt.plot(np.arange(1,7,1),B_test2,'--yD',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 2')
plt.errorbar(np.arange(1,7,1),B_test2,yerr=0.5,fmt='',linestyle="None",color='y')
plt.plot(np.arange(1,7,1),B_test1,'--mX',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 1')
plt.errorbar(np.arange(1,7,1),B_test1,yerr=0.5,fmt='',linestyle="None",color='m')
plt.ylim(0,25)
plt.xlim(0,7)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7],
           [r'$top$', r'$1$', r'$2$', r'$3$',r'$4$', r'$5$', r'$6$', r'$bottom$'])
plt.xlabel('Circuit number')
plt.ylabel('$T_{sup}$ [$^{\circ}$C]')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame  = leg.get_frame()  
frame.set_linewidth(0.5)
plt.savefig('superheat/T_sup_baseline.pdf')
plt.show()


##########plot superheat -- Modified##########
plt.plot(np.arange(1,7,1),M_testC,'--ko',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test C')
plt.errorbar(np.arange(1,7,1),M_testC,yerr=0.5,fmt='',linestyle="None",color='k')
plt.plot(np.arange(1,7,1),M_testB,'--bs',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test B')
plt.errorbar(np.arange(1,7,1),M_testB,yerr=0.5,fmt='',linestyle="None",color='b')
plt.plot(np.arange(1,7,1),M_test6,'--r^',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 6')
plt.errorbar(np.arange(1,7,1),M_test6,yerr=0.5,fmt='',linestyle="None",color='r')
plt.plot(np.arange(1,7,1),M_test5,'--g*',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 5')
plt.errorbar(np.arange(1,7,1),M_test5,yerr=0.5,fmt='',linestyle="None",color='g')
plt.plot(np.arange(1,7,1),M_test4,'--P',markersize=5,markeredgewidth=0.1,alpha=0.9,color='brown',label=r'Test 4')
plt.errorbar(np.arange(1,7,1),M_test4,yerr=0.5,fmt='',linestyle="None",color='brown')
plt.plot(np.arange(1,7,1),M_test3,'--cH',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 3')
plt.errorbar(np.arange(1,7,1),M_test3,yerr=0.5,fmt='',linestyle="None",color='c')
plt.plot(np.arange(1,7,1),M_test2,'--yD',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 2')
plt.errorbar(np.arange(1,7,1),M_test2,yerr=0.5,fmt='',linestyle="None",color='y')
plt.plot(np.arange(1,7,1),M_test1,'--mX',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 1')
plt.errorbar(np.arange(1,7,1),M_test1,yerr=0.5,fmt='',linestyle="None",color='m')
plt.ylim(0,25)
plt.xlim(0,7)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7],
           [r'$top$', r'$1$', r'$2$', r'$3$',r'$4$', r'$5$', r'$6$', r'$bottom$'])
plt.xlabel('Circuit number')
plt.ylabel('$T_{sup}$ [$^{\circ}$C]')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame  = leg.get_frame()  
frame.set_linewidth(0.5)
plt.savefig('superheat/T_sup_modified.pdf')
plt.show()  


##########plot superheat -- Interleaved##########
plt.plot(np.arange(1,7,1),I_testC,'--ko',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test C')
plt.errorbar(np.arange(1,7,1),I_testC,yerr=0.5,fmt='',linestyle="None",color='k')
plt.plot(np.arange(1,7,1),I_testB,'--bs',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test B')
plt.errorbar(np.arange(1,7,1),I_testB,yerr=0.5,fmt='',linestyle="None",color='b')
plt.plot(np.arange(1,7,1),I_test6,'--r^',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 6')
plt.errorbar(np.arange(1,7,1),I_test6,yerr=0.5,fmt='',linestyle="None",color='r')
plt.plot(np.arange(1,7,1),I_test5,'--g*',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 5')
plt.errorbar(np.arange(1,7,1),I_test5,yerr=0.5,fmt='',linestyle="None",color='g')
plt.plot(np.arange(1,7,1),I_test4,'--P',markersize=5,markeredgewidth=0.1,alpha=0.9,color='brown',label=r'Test 4')
plt.errorbar(np.arange(1,7,1),I_test4,yerr=0.5,fmt='',linestyle="None",color='brown')
plt.plot(np.arange(1,7,1),I_test3,'--cH',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 3')
plt.errorbar(np.arange(1,7,1),I_test3,yerr=0.5,fmt='',linestyle="None",color='c')
plt.plot(np.arange(1,7,1),I_test2,'--yD',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 2')
plt.errorbar(np.arange(1,7,1),I_test2,yerr=0.5,fmt='',linestyle="None",color='y')
plt.plot(np.arange(1,7,1),I_test1,'--mX',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 1')
plt.errorbar(np.arange(1,7,1),I_test1,yerr=0.5,fmt='',linestyle="None",color='m')
plt.ylim(0,25)
plt.xlim(0,7)
plt.xticks([0, 1, 2, 3, 4, 5, 6, 7],
           [r'$top$', r'$1$', r'$2$', r'$3$',r'$4$', r'$5$', r'$6$', r'$bottom$'])
plt.xlabel('Circuit number')
plt.ylabel('$T_{sup}$ [$^{\circ}$C]')
leg = plt.legend(loc='best',fancybox=False,numpoints=1)
frame  = leg.get_frame()  
frame.set_linewidth(0.5)
plt.savefig('superheat/T_sup_interleaved.pdf')
plt.show()   
