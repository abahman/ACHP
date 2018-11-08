'''
Created on Oct 26, 2017

@author: ammarbahman

Note: this file plots the Ts/Ph diagrams for gas cooler for GL 2018 paper

'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CoolProp as CP
from CoolProp.Plots import PropertyPlot
from CoolProp.CoolProp import PropsSI
import matplotlib as mpl
mpl.style.use('classic')
mpl.style.use('Elsevier.mplstyle')
mpl.rcParams['mathtext.fontset'] = 'custom'

#===============================================================================
# Latex render
#===============================================================================
# #mpl.use('pgf')
#  
# def figsize(scale):
#     fig_width_pt = 469.755                          # Get this from LaTeX using \the\textwidth
#     inches_per_pt = 1.0/72.27                       # Convert pt to inch
#     golden_mean = (np.sqrt(5.0)-1.0)/2.0            # Aesthetic ratio (you could change this)
#     fig_width = fig_width_pt*inches_per_pt*scale    # width in inches
#     fig_height = fig_width*golden_mean              # height in inches
#     fig_size = [fig_width,fig_height]
#     return fig_size
#  
# pgf_with_latex = {                      # setup matplotlib to use latex for output
# "pgf.texsystem": "pdflatex",        # change this if using xetex or lautex
# "text.usetex": True,                # use LaTeX to write all text
# "font.family": "serif",
# "mathtext.fontset" : "custom",
# "font.serif": [],                   # blank entries should cause plots to inherit fonts from the document
# "font.sans-serif": [],
# "font.monospace": [],
# "axes.labelsize": 10,               # LaTeX default is 10pt font.
# "font.size": 10,
# "legend.fontsize": 8,               # Make the legend/label fonts a little smaller
# "legend.labelspacing":0.2,
# "xtick.labelsize": 8,
# "ytick.labelsize": 8,
# "figure.figsize": figsize(0.9),     # default fig size of 0.9 textwidth
# "pgf.preamble": [
# r"\usepackage[utf8x]{inputenc}",    # use utf8 fonts becasue your computer can handle it :)
# r"\usepackage[T1]{fontenc}",        # plots will be generated using this preamble
#         ]
#     }
# mpl.rcParams.update(pgf_with_latex)
#===============================================================================
# END of Latex render
#===============================================================================

#Case 1 (supercritical -- liquid m=0.1)
T = np.array([123.1,30.9782,26.7684553578])
P = np.array([11000.0,11000.0,11000.0])
h = np.array([530.718532137,271.140002247,258.961045382])
s = np.array([1.98539276298,1.21620754604,1.1758863373])
Tair = np.array([15.0,41.1110409173])
sair = np.array([1.175886337,1.98539276298])

#Case 2 (supercritical only m=0.2)
# T = np.array([123.1,42.7627777977])
# P = np.array([11000.0,11000.0])
# h = np.array([530.718532137,314.409614245])
# s = np.array([1.98539276298,1.3556025177])
# Tair = np.array([15.0,56.5668475102])
# sair = np.array([1.3556025177,1.98539276298])



ref_fluid = 'HEOS::R744'
SS = []
TT = np.linspace(T[0],T[-1],num=1000)
for t in TT:
    SS.append(PropsSI("S","P",P[0]*1000,"T",t+273.15,ref_fluid)/1000)

Tcrit = PropsSI("Tcrit",ref_fluid) #[K]
Pcrit = PropsSI("Pcrit",ref_fluid) #[Pa]
hcrit = PropsSI("H","P",Pcrit,"T",Tcrit,ref_fluid) #[J/kg]
scrit = PropsSI("S","P",Pcrit,"T",Tcrit,ref_fluid) #[J/kg-K]

Ptriple = PropsSI("PTRIPLE",ref_fluid) #[Pa]
Ttriple = PropsSI("TTRIPLE",ref_fluid) #[K]

#Plot P-h diagram
from pylab import rcParams
rcParams['figure.figsize'] = 6.5, 5
ph_plot_R744 = PropertyPlot(ref_fluid, 'Ph',unit_system='KSI')
ph_plot_R744._set_axis_limits([50000.0, 900000.0, 50000.0, 200000000.0]) #the axis limit are changed so that the CoolProp plot function can extrapolate the lines to supercritical region
ph_plot_R744.calc_isolines(CP.iT, iso_range=[Tcrit],num=1,points=500)
ph_plot_R744.props[CP.iT]['color'] = 'black'
ph_plot_R744.props[CP.iT]['lw'] = '1.0'
ph_plot_R744.draw()
ph_plot_R744.isolines.clear()
ph_plot_R744._set_axis_limits([50000.0, 900000.0, 50000.0, 200000000.0])
ph_plot_R744.calc_isolines(CP.iQ, iso_range=[0,1],num=2,points=500)
ph_plot_R744.props[CP.iQ]['lw'] = 0.5
ph_plot_R744.draw()
ph_plot_R744.isolines.clear()
ph_plot_R744.title('P-h R744')
ph_plot_R744.xlabel(r'$h$ [kJ kg$^{-1}$]')
ph_plot_R744.ylabel(r'$P$ [kPa]')
ph_plot_R744.axis.set_yscale('log')
#ph_plot_R407C.grid()
plt.axhline(y=Pcrit/1000, color='k', linestyle='-',linewidth=1.0)
plt.plot(hcrit/1000,Pcrit/1000,'^g',markersize=8,label='Critical point')

plt.plot(h,P,'or-',markersize=8,label='Refrigerant side')
plt.title('')
plt.xlim(0,600)
plt.xticks([0,100,200,300,400,500,600],
           [r'0',r'100',r'200',r'300',r'400',r'500',r'600'])
plt.ylim(600,100000)
#plt.yticks([600,1000,10000,100000],
#           [r'600',r'1000',r'10000',r'100000'])
plt.text(550,80000,'R-744',ha='center',va='top')
plt.text(450,30000,'Supercritical',ha='center',va='top')
plt.text(470,1200,'Vapor',ha='center',va='top')
plt.text(535,5000,'Supercritical\n vapor',ha='center',va='top')
plt.text(130,30000,'Supercritical liquid',ha='center',va='top')
plt.text(80,3000,'Liquid',ha='center',va='top')
plt.text(300,2000,'Two phase',ha='center',va='top')
plt.text(40,9.5e3,'$P=P_{cr}$',ha='center',va='top')
plt.text(245,9e4,'$T=T_{cr}$',ha='center',va='top',rotation=82)
leg=plt.legend(loc='upper left',fancybox=False,numpoints=1)
frame=leg.get_frame()  
frame.set_linewidth(1.0)
plt.tight_layout() 
ph_plot_R744.savefig('Ph_case1.pdf')    
ph_plot_R744.show()
plt.close()

#Plot T-s diagram
ts_plot_R744 = PropertyPlot(ref_fluid, 'Ts',unit_system='EUR') #'EUR' is bar, kJ, C; 'KSI' is kPa, kJ, K; 'SI' is Pa, J, K
ts_plot_R744.calc_isolines(CP.iP, iso_range=[74773/1000], num=1,points=500)
ts_plot_R744.props[CP.iP]['color'] = 'black'
ts_plot_R744.props[CP.iP]['lw'] = '1.0'
ts_plot_R744.draw()
ts_plot_R744.isolines.clear()
ts_plot_R744.calc_isolines(CP.iQ, iso_range=[0,1],num=2,points=500)
ts_plot_R744.props[CP.iQ]['lw'] = 0.5
ts_plot_R744.draw()
ts_plot_R744.isolines.clear()
ts_plot_R744.title('T-s R744')
ts_plot_R744.xlabel(r'$s$ [kJ kg$^{-1}$ K$^{-1}$]')
ts_plot_R744.ylabel(r'$T$ [$\degree$C]')
#ts_plot_R407C.grid()
plt.axhline(y=Tcrit-273.15, color='k', linestyle='-',linewidth=1.0)
plt.plot(scrit/1000,Tcrit-273.15,'^g',markersize=8,label='Critical point')

plt.plot(sair,Tair,'sb-',markersize=8,label='Air side')
plt.plot(SS,TT,'or-',markersize=8,markevery=[0,-1],label='Refrigerant side')
plt.plot(s,T,'or',markersize=8)
plt.title('')
plt.ylim([-50,150])
plt.yticks([-50,-25,0,25,50,75,100,125,150],
           [r'$-50$',r'$-25$',r'0',r'25',r'50',r'75',r'100',r'125',r'150'])
plt.xlim([0.5,2.5])
plt.text(2.35,143,'R-744',ha='center',va='top')
plt.text(1.2,90,'Supercritical',ha='center',va='top')
plt.text(2.2,10,'Vapor',ha='center',va='top')
plt.text(2.225,80,'Supercritical\n vapor',ha='center',va='top')
plt.text(0.75,20,'Supercritical\n liquid',ha='center',va='top')
plt.text(1,-20,'Liquid',ha='center',va='top')
plt.plot([0.9, 1], [-10, -20],'k',linewidth=0.2)
plt.text(1.5,0,'Two phase',ha='center',va='top')
plt.text(2.09,147,'$P=P_{cr}$',ha='center',va='top',rotation=75)
plt.text(0.63,41,'$T=T_{cr}$',ha='center',va='top')
leg=plt.legend(loc='upper left',fancybox=False, numpoints=1)
frame=leg.get_frame()  
frame.set_linewidth(1.0)
plt.tight_layout() 
ts_plot_R744.savefig('Ts_case1.pdf')    
ts_plot_R744.show()
