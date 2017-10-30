'''
Created on Oct 26, 2017

@author: ammarbahman

Note: this file plots the Ts/Ph diagrams for gas cooler for GL 2018 paper

'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import CoolProp
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

# df = pd.read_excel("Ts_Ph_test1.xlsx")
# B_P = np.array(df[:12]["P"])
# B_T = np.array(df[:12]["T"])
# B_h = np.array(df[:12]["h"])
# B_s = np.array(df[:12]["s"])
# 
# M_P = np.array(df[16:28]["P"])
# M_T = np.array(df[16:28]["T"])
# M_h = np.array(df[16:28]["h"])
# M_s = np.array(df[16:28]["s"])
# 
# I_P = np.array(df[32:44]["P"])
# I_T = np.array(df[32:44]["T"])
# I_h = np.array(df[32:44]["h"])
# I_s = np.array(df[32:44]["s"])
# 
# df = pd.read_excel("x=0.8.xlsx")
# Px8 = np.array(df[:]["p"])
# hx8 = np.array(df[:]["h"])
# Tx8 = np.array(df[:]["T"])
# sx8 = np.array(df[:]["s"])
# 
# df = pd.read_excel("x=0.6.xlsx")
# Px6 = np.array(df[:]["p"])
# hx6 = np.array(df[:]["h"])
# Tx6 = np.array(df[:]["T"])
# sx6 = np.array(df[:]["s"])
# 
# df = pd.read_excel("x=0.4.xlsx")
# Px4 = np.array(df[:]["p"])
# hx4 = np.array(df[:]["h"])
# Tx4 = np.array(df[:]["T"])
# sx4 = np.array(df[:]["s"])
# 
# df = pd.read_excel("x=0.2.xlsx")
# Px2 = np.array(df[:]["p"])
# hx2 = np.array(df[:]["h"])
# Tx2 = np.array(df[:]["T"])
# sx2 = np.array(df[:]["s"])
# 
# df = pd.read_excel("x=0.xlsx")
# Tx0 = np.array(df[:]["T"])
# sx0 = np.array(df[:]["s"])
# 
# df = pd.read_excel("x=1.xlsx")
# Tx1 = np.array(df[:]["T"])
# sx1 = np.array(df[:]["s"])
               
ref_fluid = 'HEOS::R744'
Tcrit = PropsSI("Tcrit",ref_fluid) #[K]
Pcrit = PropsSI("Pcrit",ref_fluid) #[Pa]
hcrit = PropsSI("H","P",Pcrit,"T",Tcrit,ref_fluid) #[J/kg]
scrit = PropsSI("S","P",Pcrit,"T",Tcrit,ref_fluid) #[J/kg-K]

Ptriple = PropsSI("PTRIPLE",ref_fluid) #[Pa]
Ttriple = PropsSI("TTRIPLE",ref_fluid) #[K]

#Plot P-h diagram 
ph_plot_R744 = PropertyPlot(ref_fluid, 'Ph',unit_system='KSI')
ph_plot_R744._set_axis_limits([50000.0, 900000.0, 50000.0, 200000000.0]) #the axis limit are changed so that the CoolProp plot function can extrapolate the lines to supercritical region
ph_plot_R744.calc_isolines(CoolProp.iT, iso_range=[Tcrit],num=1,points=500)
ph_plot_R744.props[CoolProp.iT]['color'] = 'black'
ph_plot_R744.props[CoolProp.iT]['lw'] = '1.0'
ph_plot_R744.draw()
ph_plot_R744.isolines.clear()
ph_plot_R744._set_axis_limits([50000.0, 900000.0, 50000.0, 200000000.0])
ph_plot_R744.calc_isolines(CoolProp.iQ, iso_range=[0,1],num=2,points=500)
ph_plot_R744.title('P-h R744')
ph_plot_R744.xlabel(r'$h$ [kJ/kg]')
ph_plot_R744.ylabel(r'$P$ [kPa]')
ph_plot_R744.axis.set_yscale('log')
#ph_plot_R407C.grid()
plt.axhline(y=Pcrit/1000, color='k', linestyle='-',linewidth=1.0)
plt.plot(hcrit/1000,Pcrit/1000,'.r',markersize=10,label='Critical point')
#plt.plot(M_h,M_P,'r--',linewidth=1.5, label='Modified')
#plt.plot(I_h,I_P,'g:',linewidth=1.5,label='Interleaved')
#plt.plot(hx8,Px8,color="grey",linewidth=0.25)
#plt.text(364, 600, 'x=0.8',color="grey",fontsize=5,rotation=70)
#plt.plot(hx6,Px6,color="grey",linewidth=0.25)
#plt.text(322, 600, 'x=0.6',color="grey",fontsize=5,rotation=67)
#plt.plot(hx4,Px4,color="grey",linewidth=0.25)
#plt.text(279, 600, 'x=0.4',color="grey",fontsize=5,rotation=64)
#plt.plot(hx2,Px2,color="grey",linewidth=0.25)
plt.title('')
plt.xlim(0,600)
plt.xticks([0,100,200,300,400,500,600],
           [r'0',r'100',r'200',r'300',r'400',r'500',r'600'])
plt.ylim(600,100000)
plt.yticks([600,1000,10000,100000],
           [r'600',r'1000',r'10000',r'100000'])
plt.text(550,80000,'R-744',ha='center',va='top')
plt.text(450,30000,'Supercritical',ha='center',va='top')
plt.text(470,1200,'Vapor',ha='center',va='top')
plt.text(540,4000,'Supercritical\n vapor',ha='center',va='top')
plt.text(130,30000,'Supercritical liquid',ha='center',va='top')
plt.text(80,3000,'Liquid',ha='center',va='top')
plt.text(300,2000,'Two phase',ha='center',va='top')
leg=plt.legend(loc='upper left',fancybox=False,numpoints=1)
frame=leg.get_frame()  
frame.set_linewidth(1.0)
plt.tight_layout() 
ph_plot_R744.savefig('p-h.pdf')    
ph_plot_R744.show()
plt.close()

#Plot T-s diagram
ts_plot_R744 = PropertyPlot(ref_fluid, 'Ts',unit_system='EUR') #'EUR' is bar, kJ, C; 'KSI' is kPa, kJ, K; 'SI' is Pa, J, K
ts_plot_R744.calc_isolines(CoolProp.iP, iso_range=[74773/1000], num=1,points=500)
ts_plot_R744.props[CoolProp.iP]['color'] = 'black'
ts_plot_R744.props[CoolProp.iP]['lw'] = '1.0'
ts_plot_R744.draw()
ts_plot_R744.isolines.clear()
ts_plot_R744.calc_isolines(CoolProp.iQ, iso_range=[0,1],num=2,points=500)
ts_plot_R744.title('T-s R744')
ts_plot_R744.xlabel(r'$s$ [kJ/kg-K]')
ts_plot_R744.ylabel(r'$T$ [$\degree$C]')
#ts_plot_R407C.grid()
plt.plot(scrit/1000,Tcrit-273.15,'.r',markersize=10,label='Critical point')
plt.axhline(y=Tcrit-273.15, color='k', linestyle='-',linewidth=1.0)
#plt.plot(M_s,M_T,'r--',linewidth=1.5, label='Modified')
#plt.plot(I_s,I_T,'g:',linewidth=1.5,label='Interleaved')
#plt.plot(sx0,Tx0,color="k",linewidth=0.6)
#plt.plot(sx1,Tx1,color="k",linewidth=0.6)
#plt.plot(sx8,Tx8,color="grey",linewidth=0.25)
#plt.text(1.608, -5, 'x=0.8',color="grey",fontsize=5,rotation=90)
#plt.plot(sx6,Tx6,color="grey",linewidth=0.25)
#plt.text(1.431, -5, 'x=0.6',color="grey",fontsize=5,rotation=65)
#plt.plot(sx4,Tx4,color="grey",linewidth=0.25)
#plt.text(1.26, -5, 'x=0.4',color="grey",fontsize=5,rotation=48)
#plt.plot(sx2,Tx2,color="grey",linewidth=0.25)
plt.title('')
plt.ylim([-50,100])
plt.yticks([-50,-25,0,25,50,75,100],
           [r'-50',r'-25',r'0',r'25',r'50',r'75',r'100',r'125'])
plt.xlim([0.5,2.5])
plt.text(2.35,120,'R-744',ha='center',va='top')
plt.text(1.2,70,'Supercritical',ha='center',va='top')
plt.text(2.2,10,'Vapor',ha='center',va='top')
plt.text(2.2,65,'Supercritical\n vapor',ha='center',va='top')
plt.text(0.75,20,'Supercritical\n liquid',ha='center',va='top')
plt.text(1,-20,'Liquid',ha='center',va='top')
plt.plot([0.9, 1], [-10, -20],'k',linewidth=0.2)
plt.text(1.5,0,'Two phase',ha='center',va='top')
leg=plt.legend(loc='upper left',fancybox=False, numpoints=1)
frame=leg.get_frame()  
frame.set_linewidth(1.0)
plt.tight_layout() 
ts_plot_R744.savefig('T-s.pdf')    
ts_plot_R744.show()


    
    
##########plot superheat -- Baseline##########
# plt.plot(np.arange(1,7,1),B_testC,'--ko',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test C')
# plt.errorbar(np.arange(1,7,1),B_testC,yerr=0.5,fmt='',linestyle="None",color='k')
# plt.plot(np.arange(1,7,1),B_testB,'--bs',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test B')
# plt.errorbar(np.arange(1,7,1),B_testB,yerr=0.5,fmt='',linestyle="None",color='b')
# plt.plot(np.arange(1,7,1),B_test6,'--r^',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 6')
# plt.errorbar(np.arange(1,7,1),B_test6,yerr=0.5,fmt='',linestyle="None",color='r')
# plt.plot(np.arange(1,7,1),B_test5,'--g*',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 5')
# plt.errorbar(np.arange(1,7,1),B_test5,yerr=0.5,fmt='',linestyle="None",color='g')
# plt.plot(np.arange(1,7,1),B_test4,'--P',markersize=5,markeredgewidth=0.1,alpha=0.9,color='brown',label=r'Test 4')
# plt.errorbar(np.arange(1,7,1),B_test4,yerr=0.5,fmt='',linestyle="None",color='brown')
# plt.plot(np.arange(1,7,1),B_test3,'--cH',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 3')
# plt.errorbar(np.arange(1,7,1),B_test3,yerr=0.5,fmt='',linestyle="None",color='c')
# plt.plot(np.arange(1,7,1),B_test2,'--yD',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 2')
# plt.errorbar(np.arange(1,7,1),B_test2,yerr=0.5,fmt='',linestyle="None",color='y')
# plt.plot(np.arange(1,7,1),B_test1,'--mX',markersize=5,markeredgewidth=0.1,alpha=0.9,label=r'Test 1')
# plt.errorbar(np.arange(1,7,1),B_test1,yerr=0.5,fmt='',linestyle="None",color='m')
# plt.ylim(0,25)
# plt.xlim(0,7)
# plt.xticks([0, 1, 2, 3, 4, 5, 6, 7],
#            [r'$top$', r'$1$', r'$2$', r'$3$',r'$4$', r'$5$', r'$6$', r'$bottom$'])
# plt.xlabel('Circuit number')
# plt.ylabel('$T_{sup}$ [$^{\circ}$C]')
# leg = plt.legend(loc='best',fancybox=False,numpoints=1)
# frame  = leg.get_frame()  
# frame.set_linewidth(0.5)
# plt.savefig('T_sup_baseline.pdf')
# plt.show()


